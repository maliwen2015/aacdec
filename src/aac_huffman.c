/*
 * aacdec - AAC Huffman decoding.
 *
 * The codeword bit patterns are standardized in ISO/IEC 13818-7 and
 * ISO/IEC 14496-3.  The decoder itself builds its own binary decode trees at
 * start-up from the code lengths and codewords.
 */
#include "aac_huffman.h"

#include <stdlib.h>
#include <string.h>

#include "aac_tables_data.h"

typedef struct
{
    int16_t child[2];
    int16_t symbol;            /* >= 0 for a leaf, -1 for an internal node */
} huff_node;

typedef struct
{
    huff_node* nodes;
    int count;
    int cap;
    int symbols;
} huff_tree;

static const uint16_t* const spectral_codes[11] = {
    aac_codes_1, aac_codes_2, aac_codes_3, aac_codes_4,
    aac_codes_5, aac_codes_6, aac_codes_7, aac_codes_8,
    aac_codes_9, aac_codes_10, aac_codes_11
};

static const uint8_t* const spectral_bits[11] = {
    aac_bits_1, aac_bits_2, aac_bits_3, aac_bits_4,
    aac_bits_5, aac_bits_6, aac_bits_7, aac_bits_8,
    aac_bits_9, aac_bits_10, aac_bits_11
};

static const int spectral_sizes[11] = {
    81, 81, 81, 81, 81, 81, 64, 64, 169, 169, 289
};

static huff_tree trees[11];
static huff_tree sf_tree;
static int huffman_ready;

static int tree_add(huff_tree* t, int symbol, uint32_t code, int len)
{
    int node = 0;
    int i;

    for (i = len - 1; i >= 0; i--) {
        int bit = (int)((code >> i) & 1u);
        int next;

        if (t->nodes[node].symbol >= 0)
            return -1; /* prefix collision */

        next = t->nodes[node].child[bit];
        if (next < 0) {
            if (t->count >= t->cap)
                return -1;
            next = t->count++;
            t->nodes[next].child[0] = -1;
            t->nodes[next].child[1] = -1;
            t->nodes[next].symbol = -1;
            t->nodes[node].child[bit] = (int16_t)next;
        }
        node = next;
    }

    if (t->nodes[node].symbol >= 0 ||
        t->nodes[node].child[0] >= 0 || t->nodes[node].child[1] >= 0)
        return -1;
    t->nodes[node].symbol = (int16_t)symbol;
    return 0;
}

static int tree_build_init(huff_tree* t, int symbols)
{
    memset(t, 0, sizeof(*t));
    t->cap = symbols * 2 + 2;
    t->nodes = (huff_node*)malloc((size_t)t->cap * sizeof(*t->nodes));
    if (t->nodes == NULL)
        return -1;
    t->count = 1;
    t->nodes[0].child[0] = -1;
    t->nodes[0].child[1] = -1;
    t->nodes[0].symbol = -1;
    t->symbols = symbols;
    return 0;
}

static int tree_build16(huff_tree* t, const uint16_t* codes,
                        const uint8_t* bits, int symbols)
{
    int i;
    if (tree_build_init(t, symbols) < 0)
        return -1;
    for (i = 0; i < symbols; i++) {
        if (tree_add(t, i, codes[i], bits[i]) < 0) {
            free(t->nodes);
            t->nodes = NULL;
            return -1;
        }
    }
    return 0;
}

static int tree_build32(huff_tree* t, const uint32_t* codes,
                        const uint8_t* bits, int symbols)
{
    int i;
    if (tree_build_init(t, symbols) < 0)
        return -1;
    for (i = 0; i < symbols; i++) {
        if (tree_add(t, i, codes[i], bits[i]) < 0) {
            free(t->nodes);
            t->nodes = NULL;
            return -1;
        }
    }
    return 0;
}

static int tree_decode(aac_bits* b, const huff_tree* t)
{
    int node = 0;

    while (node >= 0 && t->nodes[node].symbol < 0) {
        int bit;
        if (aac_bits_left(b) == 0) {
            b->error = 1;
            return -1;
        }
        bit = (int)aac_bits_read(b, 1);
        node = t->nodes[node].child[bit];
    }
    if (node < 0)
        return -1;
    return t->nodes[node].symbol;
}

int aac_huffman_init(void)
{
    int i;

    if (huffman_ready)
        return 0;

    for (i = 0; i < 11; i++) {
        if (tree_build16(&trees[i], spectral_codes[i], spectral_bits[i],
                         spectral_sizes[i]) < 0)
            return -1;
    }
    if (tree_build32(&sf_tree, aac_sf_code, aac_sf_bits, 121) < 0)
        return -1;

    huffman_ready = 1;
    return 0;
}

int aac_huff_decode(aac_bits* b, int cb)
{
    if (cb < 1 || cb > 11 || !huffman_ready)
        return -1;
    return tree_decode(b, &trees[cb - 1]);
}

int aac_huff_decode_sf(aac_bits* b)
{
    if (!huffman_ready)
        return -1;
    return tree_decode(b, &sf_tree);
}
