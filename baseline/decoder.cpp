#include <decoder.h>
#include <glog/logging.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <istream>
#include <stdexcept>
#include <utility>
#include "fft.h"
#include "huffman.h"
#include <iostream>

/*
enum class Marker : uint16_t {
    SOI = 0xFFD8,
    EOI = 0xFFD9,
    COM = 0XFFFE,
    APPS = 0xFFE0,
    APPE = 0xFFEF,
    DQT = 0xFFDB,
    SOF0 = 0xFFC0,
    DHT = 0xFFC4,
    SOS = 0xFFDA
};
*/
struct ComponentInfo {
    uint16_t id;
    uint16_t H;
    uint16_t V;
    uint16_t quant_id;
};

struct ScanComponent {
    uint8_t cid;
    uint8_t dc_id;
    uint8_t ac_id;
};

class BitReader {
public:
    explicit BitReader(std::istream& in) : in_(in) {
    }

    uint16_t ReadBits(int n) {
        if (n < 1 || n > 16) {
            throw std::runtime_error("from 1 to 16 pls");
        }

        while (bits_in_buf_ < n) {
            Fillbuf();
        }

        uint16_t res = static_cast<uint16_t>(buf_ >> (bits_in_buf_ - n));
        bits_in_buf_ -= n;
        uint32_t mask = (1u << bits_in_buf_) - 1u;
        buf_ &= mask;
        return res;
    }

    int16_t ReadSigned(int n) {
        if (n == 0) {
            return 0;
        }

        uint16_t raw = ReadBits(n);
        int16_t lim = 1 << (n - 1);
        if (raw < lim) {
            raw = static_cast<uint16_t>(raw - ((1 << n) - 1));
        }
        return static_cast<int16_t>(raw);
    }

    bool Eof() const {
        return eof_;
    }

private:
    void Fillbuf() {
        if (eof_) {
            throw std::runtime_error("End of input not expected");
        }

        int nextb = Getnext();
        buf_ = (buf_ << 8) | (nextb & 0xFF);
        bits_in_buf_ += 8;
    }

    int Getnext() {
        int val = in_.get();
        if (val == EOF) {
            eof_ = true;
            throw std::runtime_error("End?");
        }

        if ((val & 0xFF) == 0xFF) {
            int nb = in_.peek();
            if (nb == EOF) {
                eof_ = true;
                throw std::runtime_error("EOF after 0xFF xd");
            }

            if ((nb & 0xFF) == 0x00) {
                in_.ignore(1);
                return 0xFF;
            } else {
                eof_ = true;
                throw std::runtime_error("Found marker");
            }
        }
        return val & 0xFF;
    }

    std::istream& in_;
    uint32_t buf_ = 0;
    int bits_in_buf_ = 0;
    bool eof_ = false;
};

/*

const std::vector <int> order =
    {0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32,
    25, 18, 11, 4, 5, 12, 19, 26, 33,
     40, 48, 41, 34, 27, 20, 13, 6, 7,
      14, 21, 28, 35, 42,  49, 56, 57, 50,
      43, 36, 29, 22, 15, 23, 30, 37, 44,
      51, 58, 59, 52, 45, 38, 31, 39, 46,
      53, 60, 61, 54, 47, 55, 62, 63};
*/

void ZigZag(std::vector<int>& bef, std::vector<int>& res) {
    std::vector<int> order = {0,  1,  8,  16, 9,  2,  3,  10, 17, 24, 32, 25, 18, 11, 4,  5,
                              12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6,  7,  14, 21, 28,
                              35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
                              58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63};

    if (bef.size() != 64) {
        throw std::runtime_error("NOT 8x8 data");
    }

    for (int i = 0; i < 64; ++i) {
        res[order[i]] = bef[i];
    }
}

static uint16_t Readu16be(std::istream& in) {
    if (in.eof()) {
        throw std::runtime_error("EOF??");
    }
    uint8_t hi = static_cast<uint8_t>(in.get());
    if (in.eof()) {
        throw std::runtime_error("EOF??");
    }
    uint8_t lo = static_cast<uint8_t>(in.get());
    return static_cast<uint16_t>((hi << 8) | lo);
}

static uint8_t Read8(std::istream& in) {
    if (in.eof()) {
        throw std::runtime_error("EOF??");
    }
    return static_cast<uint8_t>(in.get());
}

void Decodeblock(BitReader& br, HuffmanTree& dctree, HuffmanTree& actree, int& prevdc,
                 const std::vector<int>& quant, std::vector<int>& res) {
    int dclen = 0;
    while (true) {
        uint16_t b = br.ReadBits(1);
        bool mv = (b == 1);
        int val = 0;
        if (dctree.Move(mv, val)) {
            dclen = val;
            break;
        }
    }

    int16_t dcdif = br.ReadSigned(dclen);
    int dccur = prevdc + dcdif;
    prevdc = dccur;
    std::vector<int> bl(64, 0);
    bl[0] = dccur;

    for (int i = 1; i < 64;) {
        int rval = 0;
        while (true) {
            uint16_t b = br.ReadBits(1);
            bool mv = (b == 1);
            int val = 0;
            if (actree.Move(mv, val)) {
                rval = val;
                break;
            }
        }

        int run = rval / 16;
        int sz = rval % 16;
        if (run == 0 && sz == 0) {
            break;
        }

        if (i + run > 63) {
            throw std::runtime_error("Wrong amount");
        }

        i += run;

        int16_t acval = br.ReadSigned(sz);

        bl[i] = acval;
        i++;
    }
    auto bl1 = bl;
    ZigZag(bl, bl1);
    for (int i = 0; i < 64; ++i) {
        bl[i] = bl1[i] * quant[i];
    }

    std::vector<double> indct(64), outdct(64);
    for (int i = 0; i < 64; ++i) {
        indct[i] = static_cast<double>(bl[i]);
    }

    DctCalculator dctcc(8, &indct, &outdct);
    dctcc.Inverse();

    for (int i = 0; i < 64; ++i) {
        res[i] = static_cast<int>(std::round(outdct[i]));
    }
}

Image Decode(std::istream& input) {

    Image resimg;
    uint16_t marker = Readu16be(input);
    if (marker != 0xFFD8) {
        throw std::runtime_error("NO SOI MARKER!");  // SOEVIY ******
    }

    bool eoi_found = false;
    bool was_com = false;
    bool was_sof = false;
    uint16_t height = 0;
    uint16_t width = 0;
    uint16_t n_components = 0;
    std::vector<ComponentInfo> components;
    std::vector<std::vector<int>> quanttable(4, std::vector<int>(64));

    HuffmanTree dc_tables[4];
    HuffmanTree ac_tables[4];

    std::vector<ScanComponent> scan_components;

    while (!eoi_found && !input.eof()) {
        marker = Readu16be(input);

        if (marker == 0xFFD9) {  // EOI
            eoi_found = true;
        } else if (marker == 0xFFFE) {  // COM
            if (!was_com) {
                was_com = true;
            } else {
                throw std::runtime_error("two comments");
            }
            auto len = Readu16be(input);
            if (len < 2) {
                throw std::runtime_error("COM len is too low");
            }
            len -= 2;
            std::string com;
            while (len--) {
                char c = static_cast<char>(Read8(input));
                com.push_back(c);
            }
            // std::cout << "COM: " << com << "\n";
            resimg.SetComment(com);

            // OK
        } else if ((marker & 0xFFE0) == 0xFFE0) {  // APPn
            auto len = Readu16be(input);
            if (len < 2) {
                throw std::runtime_error("APPn len is too low");
            }
            len -= 2;
            while (len--) {
                Read8(input);
            }

            // OK
        } else if (marker == 0xFFDB) {  // DQT
            auto len = Readu16be(input);
            if (len < 2) {
                throw std::runtime_error("DQT len is too low");
            }
            len -= 2;

            while (len > 0) {
                auto idf = Read8(input);
                int tp = (idf / 16);
                int id = idf % 16;
                len--;
                if (id > 3) {
                    throw std::runtime_error("WTH with DQT id??");
                }

                std::vector<int> resbef(64);
                if (tp == 0) {
                    if (len < 64) {
                        throw std::runtime_error("len xdd");
                    }
                    for (int i = 0; i < 64; ++i) {
                        auto cur = Read8(input);
                        resbef[i] = cur;
                    }
                    len -= 64;
                } else {
                    if (len < 128) {
                        throw std::runtime_error("len xdd");
                    }
                    for (int i = 0; i < 64; ++i) {
                        auto cur = Readu16be(input);
                        resbef[i] = static_cast<int>(cur);
                    }
                    len -= 128;
                }

                std::vector<int> res(64);
                ZigZag(resbef, res);
                quanttable[id] = res;
                // std::cout << "QUANT:\n";
                for (int i = 0; i < 8; ++i) {
                    for (int j = 0; j < 8; ++j) {
                        // std::cout << res[i * 8 + j] << " ";
                    }
                    // std::cout << "\n";
                }

                // SEEMS OK
            }
        } else if (marker == 0xFFC0) {  // SOF0
            if (!was_sof) {
                was_sof = true;
            } else {
                throw std::runtime_error("two SOF");
            }
            auto len = Readu16be(input);
            // std::cout << len << "\n";
            if (len < 2) {
                throw std::runtime_error("SOF0 len is too low");
            }
            len -= 2;
            if (len < 6) {
                throw std::runtime_error("SOF0 len is too low");
            }
            Read8(input);
            height = Readu16be(input);
            width = Readu16be(input);
            if (height * width <= 0) {
                throw std::runtime_error("h*w <= 0");
            }
            // std::cout << n_components << "\n";
            n_components = Read8(input);
            // std::cout << n_components << "\n";
            resimg.SetSize(width, height);
            len -= 6;
            components.clear();
            if (len < 3 * n_components) {
                throw std::runtime_error("SOF0 len is too low");
            }

            for (int i = 0; i < n_components; ++i) {
                ComponentInfo ci;
                ci.id = Read8(input);
                uint8_t hv = Read8(input);
                ci.H = hv / 16;
                ci.V = hv % 16;
                ci.quant_id = Read8(input);
                components.push_back(ci);
                len -= 3;
            }

            while (len > 0) {
                Read8(input);
                len--;
            }

            // OK

        } else if (marker == 0xFFC4) {  // DHT
            auto len = Readu16be(input);
            if (len < 2) {
                throw std::runtime_error("DHT len is too low");
            }
            len -= 2;
            while (len > 0) {
                uint8_t ht_info = Read8(input);
                len--;
                int ht_class = ht_info / 16;
                int ht_id = ht_info % 16;
                if (ht_class > 1) {
                    throw std::runtime_error("wrong class for Huffman talbe");
                }
                if (ht_id > 3) {
                    throw std::runtime_error("wrong id for Huffman talbe");
                }

                if (len < 16) {
                    throw std::runtime_error("DHT len error");
                }
                std::vector<uint8_t> code_lengths(16);
                for (int i = 0; i < 16; ++i) {
                    code_lengths[i] = Read8(input);
                }
                len -= 16;
                int ttlsymb = 0;
                for (int i = 0; i < 16; ++i) {
                    ttlsymb += code_lengths[i];
                }

                if (len < static_cast<uint16_t>(ttlsymb)) {
                    throw std::runtime_error("DHT len error");
                }

                std::vector<uint8_t> symbols(ttlsymb);
                for (int i = 0; i < ttlsymb; ++i) {
                    symbols[i] = Read8(input);
                }

                len -= ttlsymb;

                HuffmanTree tree;
                tree.Build(code_lengths, symbols);

                if (ht_class == 0) {  // DC
                    dc_tables[ht_id] = std::move(tree);
                } else {  // AC
                    ac_tables[ht_id] = std::move(tree);
                }  // AC/DC xd

                // OK
            }

        } else if (marker == 0xFFDA) {  // SOS
            auto len = Readu16be(input);

            if (len < 2) {
                throw std::runtime_error("SOS len error");
            }
            len -= 2;

            if (len < 1) {
                throw std::runtime_error("SOS len error");
            }
            uint8_t ns = Read8(input);
            len--;
            if (ns != n_components) {
                throw std::runtime_error("ns != n_comp");
            }
            scan_components.clear();
            if (len < ns * 2) {
                throw std::runtime_error("SOS len error");
            }
            for (int i = 0; i < ns; ++i) {
                ScanComponent sc;
                sc.cid = Read8(input);
                if (sc.cid > n_components) {
                    throw std::runtime_error("wrong component num");
                }
                auto acdc = Read8(input);
                sc.dc_id = acdc / 16;
                sc.ac_id = acdc % 16;
                scan_components.push_back(sc);
                if (sc.dc_id > 3 || sc.ac_id > 3) {
                    throw std::runtime_error("wrond ac/dc id");
                }
            }
            len -= ns * 2;

            // OK

            if (len < 3) {
                throw std::runtime_error("len error");
            }
            Read8(input);
            int mxq = Read8(input);
            if (mxq != 63) {
                throw std::runtime_error("not till 63");
            }
            Read8(input);

            // NACHINAETSYA MYASO xd

            if (n_components == 3) {
                // 'RGB' MODE

                std::pair<int, int> hv;

                hv = {std::max(std::max(components[0].H, components[1].H), components[2].H),
                      std::max(std::max(components[0].V, components[1].V), components[2].V)};

                BitReader br(input);
                int maxh = hv.first;
                int maxv = hv.second;
                int mcux = (width + (8 * maxh - 1)) / (8 * maxh);
                int mcuy = (height + (8 * maxv - 1)) / (8 * maxv);

                std::vector<int> ply(width * height), plcb = ply, plcr = ply;

                int prevdc[3] = {0, 0, 0};

                for (int my = 0; my < mcuy; ++my) {
                    for (int mx = 0; mx < mcux; ++mx) {
                        for (int ci = 0; ci < 3; ++ci) {
                            int hi = components[ci].H;
                            int vi = components[ci].V;
                            int qid = components[ci].quant_id;
                            int dcid = scan_components[ci].dc_id;
                            int acid = scan_components[ci].ac_id;

                            HuffmanTree& dctree = dc_tables[dcid];
                            HuffmanTree& actree = ac_tables[acid];
                            std::vector<int>& quant = quanttable[qid];

                            for (int sy = 0; sy < vi; ++sy) {
                                for (int sx = 0; sx < hi; ++sx) {
                                    std::vector<int> resqq(64, 0);
                                    Decodeblock(br, dctree, actree, prevdc[ci], quant, resqq);

                                    int gx = mx * 8 * hi + sx * 8;
                                    int gy = my * 8 * vi + sy * 8;

                                    for (int qy = 0; qy < 8; ++qy) {
                                        for (int qx = 0; qx < 8; ++qx) {
                                            int yc = gy + qy;
                                            int xc = gx + qx;
                                            if (yc < height && xc < width) {
                                                int idx = yc * width + xc;
                                                int val = resqq[qy * 8 + qx];

                                                if (ci == 0) {
                                                    ply[idx] = val;
                                                } else if (ci == 1) {
                                                    plcb[idx] = val;
                                                } else {
                                                    plcr[idx] = val;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                for (int i = 0; i < width * height; ++i) {
                    int nw = ply[i] + 128;
                    nw = std::min(nw, 255);
                    nw = std::max(nw, 0);
                    ply[i] = nw;

                    nw = plcb[i] + 128;
                    nw = std::min(nw, 255);
                    nw = std::max(nw, 0);
                    plcb[i] = nw;

                    nw = plcr[i] + 128;
                    nw = std::min(nw, 255);
                    nw = std::max(nw, 0);
                    plcr[i] = nw;
                }

                std::vector<int> plcb1(width * height), plcr1 = plcb1;

                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        int ind = y * width + x;
                        int y11 = y / maxv;
                        int x11 = x / maxh;
                        int ind11 = y11 * width + x11;
                        double yy = ply[ind];
                        double cb = plcb[ind11] - 128.0;
                        double cr = plcr[ind11] - 128.0;

                        double r = yy + (1.402 * cr);
                        double g = yy - (0.34414 * cb) - (0.71414 * cr);
                        double b = yy + (1.772 * cb);

                        int r1 = static_cast<int>(std::round(r));
                        int g1 = static_cast<int>(std::round(g));
                        int b1 = static_cast<int>(std::round(b));

                        r1 = std::min(r1, 255);
                        r1 = std::max(r1, 0);
                        g1 = std::min(g1, 255);
                        g1 = std::max(g1, 0);
                        b1 = std::min(b1, 255);
                        b1 = std::max(b1, 0);

                        resimg.SetPixel(y, x, RGB{r1, g1, b1});
                    }
                }
            } else if (n_components == 1) {
                // GRAY MODE

                int mxh = components[0].H;
                int mxv = components[0].V;
                int mcux = (width + 8 * mxh - 1) / (8 * mxh);
                int mcuy = (height + 8 * mxv - 1) / (8 * mxv);

                std::vector<int> plg(width * height, 0);
                int prevdc = 0;

                BitReader br(input);

                for (int my = 0; my < mcuy; ++my) {
                    for (int mx = 0; mx < mcux; ++mx) {
                        std::vector<int> bl(64, 0);

                        int qid = components[0].quant_id;
                        int dcid = scan_components[0].dc_id;
                        int acid = scan_components[0].ac_id;

                        Decodeblock(br, dc_tables[dcid], ac_tables[acid], prevdc, quanttable[qid],
                                    bl);

                        int gx = mx * 8;
                        int gy = my * 8;

                        for (int iy = 0; iy < 8; ++iy) {
                            for (int ix = 0; ix < 8; ++ix) {
                                int yc = gy + iy;
                                int xc = gx + ix;

                                if (yc < height && xc < width) {
                                    int ind = yc * width + xc;
                                    plg[ind] = bl[iy * 8 + ix];
                                }
                            }
                        }
                    }
                }

                for (int i = 0; i < width * height; ++i) {
                    int val = plg[i] + 128;
                    val = std::max(0, val);
                    val = std::min(255, val);
                    plg[i] = val;
                }

                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        int g = plg[y * width + x];
                        resimg.SetPixel(y, x, RGB{g, g, g});
                    }
                }
            } else {
                throw std::runtime_error("too early!");
            }
        } else {
            throw std::runtime_error("wrong marker");
        }
    }

    if (!eoi_found) {
        throw std::runtime_error("NO EOI");
    }

    return resimg;
}
