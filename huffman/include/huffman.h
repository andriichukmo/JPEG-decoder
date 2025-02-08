// Don't change this file, it is not sent to the server.

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <memory>

// HuffmanTree decoder for DHT section.
class HuffmanTree {
public:
    HuffmanTree();

    HuffmanTree(const HuffmanTree&) = delete;
    HuffmanTree& operator=(const HuffmanTree&) = delete;

    HuffmanTree(HuffmanTree&&);
    HuffmanTree& operator=(HuffmanTree&&);

    // code_lengths is the array of size no more than 16 with number of
    // terminated nodes in the Huffman tree.
    // values are the values of the terminated nodes in the consecutive
    // level order.
    void Build(const std::vector<uint8_t>& code_lengths, const std::vector<uint8_t>& values);

    // Moves the state of the Huffman tree by |bit|. If the node is terminal, returns true,
    // overwrites |value| and resets the tree state. If the node is intermediate, false is
    // returned without changing the value.
    bool Move(bool bit, int& value);

    ~HuffmanTree();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
