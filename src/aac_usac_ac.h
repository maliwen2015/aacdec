/*
 * aacdec - USAC (ISO/IEC 23003-3) context-based arithmetic decoder for
 * scalefactorband_data().
 *
 * The arithmetic decoder follows the adaptive binary arithmetic decoding
 * process of ISO/IEC 23003-3: a 16 bit interval [low, high] is maintained and
 * CDFs scaled to 2^14 drive the symbol decisions.  Quantized values are
 * decoded as 2-tuples (a, b): the most significant bits come from the MSB
 * CDFs selected through the context/pk mapping, the least significant bits
 * from one of three LSB CDFs.  The context memory is carried across frames.
 */
#ifndef AAC_USAC_AC_H
#define AAC_USAC_AC_H

#include <stddef.h>
#include <stdint.h>

#include "aac_bits.h"

#define AAC_USAC_AC_ESCAPE 16

/* Persistent context memory for one arithmetic-decoded spectrum.  The
 * 2-tuple context of the previous frame is kept so that context mapping can
 * be continued across frames (arith_reset_flag restarts it). */
typedef struct
{
    uint8_t last[513];          /* previous frame context, 2-tuple index */
    int last_len;               /* number of valid entries in last[] */
} aac_usac_ac_state;

/* Decode one arithmetically coded spectrum.  coef receives lg_max float
 * coefficients (signed quantized values, the caller applies scale factors).
 * lg coefficients are decoded; the remaining lg_max - lg are zeroed.  reset
 * must be 1 for the first window after an independent frame.  Returns 0 on
 * success, -1 on bitstream errors. */
int aac_usac_ac_decode(aac_usac_ac_state* st, aac_bits* b, float* coef,
                       int lg, int lg_max, int reset);

#endif /* AAC_USAC_AC_H */
