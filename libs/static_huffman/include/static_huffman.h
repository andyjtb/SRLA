#ifndef STATICHUFFMAN_H_INCLUDED
#define STATICHUFFMAN_H_INCLUDED

#include <stdint.h>
#include "bit_stream.h"

/* Maximum number of symbols to encode */
#define STATICHUFFMAN_MAX_NUM_SYMBOLS 256

/* Huffman tree */
struct StaticHuffmanTree {
    uint32_t num_symbols;                       /* Number of encoded symbols */
    uint32_t root_node;                         /* Index of the root node */
    struct {
        uint32_t node_0;                        /* index of left child */
        uint32_t node_1;                        /* index of right child */
    } nodes[2 * STATICHUFFMAN_MAX_NUM_SYMBOLS]; /* Tree nodes */
};

/* Huffman code */
struct StaticHuffmanCodes {
    uint32_t num_symbols;                    /* Number of encoded symbols */
    struct {
        uint32_t code;                       /* Assigned code (maximum 32 bits) */
        uint8_t bit_count;                   /* Code length */
    } codes[STATICHUFFMAN_MAX_NUM_SYMBOLS];  /* Sign of each symbol */
};

#ifdef __cplusplus
extern "C" {
#endif

/* Build the Huffman tree */
void StaticHuffman_BuildHuffmanTree(
        const uint32_t *symbol_counts, uint32_t num_symbols, struct StaticHuffmanTree *tree);

/* Create code table */
void StaticHuffman_ConvertTreeToCodes(
        const struct StaticHuffmanTree *tree, struct StaticHuffmanCodes *codes);

/* Output Huffman code */
void StaticHuffman_PutCode(
        const struct StaticHuffmanCodes *codes, struct BitStream *stream, uint32_t val);

/* Get Huffman code */
uint32_t StaticHuffman_GetCode(
        const struct StaticHuffmanTree *tree, struct BitStream *stream);

#ifdef __cplusplus
}
#endif

#endif /* STATICHUFFMAN_H_INCLUDED */
