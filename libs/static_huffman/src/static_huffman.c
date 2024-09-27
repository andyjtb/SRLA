#include "static_huffman.h"
#include <stdint.h>
#include <string.h>

/* Count normalization */
static void StaticHuffman_NormalizeSymbolCounts(
    const uint32_t *symbol_counts, uint32_t num_symbols,
    uint32_t *normalized_counts, uint32_t normalized_counts_size)
{
    uint32_t node;

    assert((symbol_counts != NULL) && (normalized_counts != NULL));
    assert((2 * num_symbols) <= normalized_counts_size);

    /* Copy frequencies to working area */
    memset(normalized_counts, 0, sizeof(uint32_t) * normalized_counts_size);
    memcpy(normalized_counts, symbol_counts, sizeof(uint32_t) * num_symbols);

    /* Avoid 0 (invalid value) */
    for (node = 0; node < num_symbols; node++) {
        if (normalized_counts[node] == 0) {
            normalized_counts[node] += 1;
        }
    }
}

/* Construct a Huffman code */
void StaticHuffman_BuildHuffmanTree(
    const uint32_t *symbol_counts, uint32_t num_symbols, struct StaticHuffmanTree *tree)
{
#define SENTINEL_NODE (2 * STATICHUFFMAN_MAX_NUM_SYMBOLS)
    uint32_t min1, min2;  /* min1 is the minimum frequency, min2 is the second minimum frequency */
    uint32_t free_node, node;
    uint32_t counts_work[(2 * STATICHUFFMAN_MAX_NUM_SYMBOLS) + 1]; /* Symbol frequency (plus one at sentinel node) */

    assert((symbol_counts != NULL) && (tree != NULL));
    assert(num_symbols > 0);
    assert(num_symbols <= STATICHUFFMAN_MAX_NUM_SYMBOLS);

    /* Set the number of symbols */
    tree->num_symbols = num_symbols;

    /* Normalize frequency counts */
    StaticHuffman_NormalizeSymbolCounts(
        symbol_counts, num_symbols, counts_work, 2 * STATICHUFFMAN_MAX_NUM_SYMBOLS);
    /* Set maximum frequency for sentinel nodes */
    counts_work[SENTINEL_NODE] = UINT32_MAX;

    /* Create parent node: Use nodes after num_symbols */
    for (free_node = num_symbols; ; free_node++) {
        /* Set index to sentinel */
        min1 = min2 = SENTINEL_NODE;

        /* Find the index that gives the 1st and 2nd smallest value */
        /*
/* Start with all nodes first, and from the next time onwards, find the index that gives the minimum value including the two nodes
* and the parent nodes */
*/
        for (node = 0; node < free_node; node++) {
            /* Only refer to the frequency of the node in question if it is not an invalid value (0) */
            if (counts_work[node] > 0) {
                if (counts_work[node] < counts_work[min1]) {
                    min2 = min1;
                    min1 = node;
                } else if (counts_work[node] < counts_work[min2]) {
                    min2 = node;
                }
            }
        }
        assert(min1 != SENTINEL_NODE);

        /* Second minimum not found */
        /* -> The nodes are grouped together until there is only one. The tree root has been found. */
        if (min2 == SENTINEL_NODE) {
            break;
        }

        /* Fill parent node with information */
        /* The frequency of a parent node is the sum of its children */
        counts_work[free_node] = counts_work[min1] + counts_work[min2];
        /* Child node frequency is invalid */
        counts_work[min1] = counts_work[min2] = 0;
        /* Record the child node index */
        tree->nodes[free_node].node_0 = min1;
        tree->nodes[free_node].node_1 = min2;
    }

    assert(free_node <= (2 * STATICHUFFMAN_MAX_NUM_SYMBOLS));

    /* After incrementing with the for statement, subtract 1 and record the root */
    tree->root_node = free_node - 1;

#undef SENTINEL_NODE
}

/* Construct a code from the Huffman tree */
static void StaticHuffman_ConvertTreeToCodesCore(
    const struct StaticHuffmanTree *tree, struct StaticHuffmanCodes *codes,
    uint32_t node, uint32_t code, uint8_t bit_count)
{
    assert(tree != NULL);
    assert(codes != NULL);

    /* The referenced node index has reached a leaf */
    if (node < tree->num_symbols) {
        /* Assign code and bit count */
        codes->codes[node].code = code;
        codes->codes[node].bit_count = bit_count;
        return;
    }

    /* Lengthen the code by 1 bit */
    code <<= 1;
    bit_count++;

    /* Follow the left leaf. The least significant bit of the code is padded with a 0 */
    StaticHuffman_ConvertTreeToCodesCore(tree, codes, tree->nodes[node].node_0, code | 0, bit_count);
    /* Follow the right leaf. Add a 1 to the least significant bit of the code. */
    StaticHuffman_ConvertTreeToCodesCore(tree, codes, tree->nodes[node].node_1, code | 1, bit_count);
}

/* Create code table */
void StaticHuffman_ConvertTreeToCodes(
    const struct StaticHuffmanTree *tree, struct StaticHuffmanCodes *codes)
{
    assert((tree != NULL) && (codes != NULL));

    /* Record the number of symbols */
    codes->num_symbols = tree->num_symbols;

    /* Start recursion from the root */
    StaticHuffman_ConvertTreeToCodesCore(tree, codes, tree->root_node, 0, 0);
}

/* Output Huffman code */
void StaticHuffman_PutCode(
    const struct StaticHuffmanCodes *codes, struct BitStream *stream, uint32_t val)
{
    assert(codes != NULL);
    assert(stream != NULL);
    assert(val < codes->num_symbols);

    BitWriter_PutBits(stream, codes->codes[val].code, codes->codes[val].bit_count);
}

/* Get Huffman code */
uint32_t StaticHuffman_GetCode(
    const struct StaticHuffmanTree *tree, struct BitStream *stream)
{
    uint32_t node, bit;

    assert(tree != NULL);
    assert(stream != NULL);

    /* Set the node to the root */
    node = tree->root_node;

    /* Traverse the tree until a leaf node is reached */
    do {
        BitReader_GetBits(stream, &bit, 1);
        node = (bit == 0) ? tree->nodes[node].node_0 : tree->nodes[node].node_1;
    } while (node >= tree->num_symbols);

    assert(node < tree->num_symbols);

    return node;
}
