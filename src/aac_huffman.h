/*
 * aacdec - AAC Huffman decoding.
 */
#ifndef AAC_HUFFMAN_H
#define AAC_HUFFMAN_H

#include "aac_bits.h"

/* Decode one codeword from spectral codebook cb (1..11).  Returns the
 * codebook symbol index, or -1 on a bitstream error. */
int aac_huff_decode(aac_bits* b, int cb);

/* Decode one codeword from the scalefactor codebook (0..120). */
int aac_huff_decode_sf(aac_bits* b);

int aac_huffman_init(void);

#endif /* AAC_HUFFMAN_H */
