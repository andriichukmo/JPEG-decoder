#include <huffman.h>
// NO MORE "IF" ANY TESTS xd
struct Node {
    bool is_leaf = false;
    int val = 0;
    Node* left = nullptr;
    Node* right = nullptr;
};

class HuffmanTree::Impl {
public:
    Node* root = nullptr;
    Node* cur = nullptr;
    bool built = false;

    Impl() {
    }

    ~Impl() {
        Freenode(root);
    }

    void Freenode(Node* node) {
        if (!node) {
            return;
        }
        Freenode(node->left);
        Freenode(node->right);
        delete node;
    }

    void Reset() {
        Freenode(root);
        root = nullptr;
        cur = nullptr;
        built = false;
    }

    void Build(const std::vector<uint8_t>& code_lengths, const std::vector<uint8_t>& values) {
        Reset();
        if (code_lengths.size() > 16) {
            throw std::invalid_argument("too large");
        }
        int max_len = 0;
        int ttl = 0;
        for (size_t i = 0; i < code_lengths.size(); ++i) {
            if (code_lengths[i] > 0) {
                max_len = static_cast<int>(i) + 1;
            }
            ttl += code_lengths[i];
        }

        if (ttl != static_cast<int>(values.size())) {
            throw std::invalid_argument("Wrong match");
        }

        {
            uint32_t code = 0;
            int prev = 0;
            for (int len = 1; len <= max_len; ++len) {
                code <<= (len - prev);
                prev = len;
                if (code + code_lengths[len - 1] > (1u << len)) {
                    throw std::invalid_argument("Wrond distribution");
                }
                code += code_lengths[len - 1];
            }
        }

        root = new Node();
        cur = root;
        int ind = 0;
        uint32_t code = 0;
        int prev = 0;
        for (int len = 1; len <= max_len; ++len) {
            code <<= (len - prev);
            prev = len;
            int cnt = code_lengths[len - 1];
            for (int i = 0; i < cnt; ++i) {
                uint8_t val = values[ind++];
                InsertCode(code, len, val);
                code++;
            }
        }
        built = true;
    }

    void InsertCode(uint32_t code, int len, uint8_t val) {
        Node* cur = root;
        for (int i = len - 1; i >= 0; --i) {
            bool bit = (code >> i) & 1;
            if (!bit) {
                if (!cur->left) {
                    cur->left = new Node();
                }
                cur = cur->left;
            } else {
                if (!cur->right) {
                    cur->right = new Node();
                }
                cur = cur->right;
            }
        }
        if (cur->is_leaf) {
            throw std::invalid_argument("duplicate");
        }
        cur->is_leaf = true;
        cur->val = val;
    }

    bool Move(bool bit, int& val) {
        if (!built) {
            throw std::invalid_argument("Not built");
        }
        if (!root) {
            throw std::invalid_argument("Empty tree");
        }
        if (!cur) {
            cur = root;
        }
        Node* nxt = (bit ? cur->right : cur->left);
        if (!nxt) {
            throw std::invalid_argument("Invalid path");
        }

        cur = nxt;
        if (cur->is_leaf) {
            val = cur->val;
            cur = root;
            return true;
        }
        return false;
    }
};

HuffmanTree::HuffmanTree() : impl_(new Impl()){};
HuffmanTree::HuffmanTree(HuffmanTree&& other) = default;
HuffmanTree& HuffmanTree::operator=(HuffmanTree&& other) = default;

HuffmanTree::~HuffmanTree() = default;

void HuffmanTree::Build(const std::vector<uint8_t>& code_lengths,
                        const std::vector<uint8_t>& values) {
    impl_->Build(code_lengths, values);
}

bool HuffmanTree::Move(bool bit, int& value) {
    return impl_->Move(bit, value);
}
