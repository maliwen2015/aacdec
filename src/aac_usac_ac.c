/*
 * aacdec - USAC (ISO/IEC 23003-3) context-based arithmetic decoder for
 * scalefactorband_data().
 */
#include <string.h>
#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#endif

#include "aac_usac_ac.h"
#include "aac_usac_ac_tables.h"
#include "aac_tables.h"

/* Arithmetic decoder interval state.  low/high/val are 16 bit values; the
 * CDFs are scaled to 2^14 (matching the normative statBits = 14). */
typedef struct
{
    uint32_t low;
    uint32_t high;
    uint32_t val;
} aac_ac_engine;

static void aac_ac_engine_init(aac_ac_engine* e, aac_bits* b)
{
    e->low = 0;
    e->high = 0xFFFFu;
    e->val = aac_bits_read(b, 16);
}

/* Decode one symbol using the CDF cdf[0..cdf_len-1] (descending cumulative
 * frequencies, last entry 0).  Returns the symbol index. */
static int aac_ac_decode(aac_ac_engine* e, aac_bits* b,
                         const uint16_t* cdf, int cdf_len)
{
    int val = (int)e->val;
    int low = (int)e->low;
    int high = (int)e->high;
    int rng = high - low + 1;
    int c = ((val - low + 1) << 14) - 1;
    int sym = 0;

    while (sym < cdf_len - 1 && cdf[sym] * rng > c)
        sym++;

    if (sym)
        high = low + ((rng * cdf[sym - 1]) >> 14) - 1;
    low += (rng * cdf[sym]) >> 14;

    /* Renormalisation.  The interval registers must stay wider than 16 bits
     * while the loop runs; truncating to 16 bits here desynchronises val from
     * low/high and corrupts long symbol runs. */
    for (;;) {
        if (high < 0x8000) {
            /* nothing */
        } else if (low >= 0x8000) {
            val -= 0x8000;
            low -= 0x8000;
            high -= 0x8000;
        } else if (low >= 0x4000 && high < 0xC000) {
            val -= 0x4000;
            low -= 0x4000;
            high -= 0x4000;
        } else {
            break;
        }
        low += low;
        high += high + 1;
        val = (val << 1) | aac_bits_read(b, 1);
    }

    e->low = (uint16_t)low;
    e->high = (uint16_t)high;
    e->val = (uint16_t)val;
    return sym;
}

/* Map the 2-tuple context to a probability model index ("pk").  The mapping
 * tables implement the normative context model; entries pair a context
 * threshold with the pk index. */
static int aac_ac_get_pk(uint32_t c)
{
    int i_min = -1;
    int i_max = (int)(sizeof(aac_usac_ac_lookup_m) /
                      sizeof(aac_usac_ac_lookup_m[0])) - 1;

    while ((i_max - i_min) > 1) {
        int i = i_min + ((i_max - i_min) / 2);
        uint32_t j = aac_usac_ac_hash_m[i];
        if (c < (j >> 8))
            i_max = i;
        else if (c > (j >> 8))
            i_min = i;
        else
            return (int)(j & 0xFF);
    }
    return aac_usac_ac_lookup_m[i_max];
}

/* arith_map_context: prepare the context memory for a new frame. */
static uint32_t aac_ac_map_process(aac_usac_ac_state* st, int reset, int N)
{
    if (reset) {
        memset(st->last, 0, sizeof(st->last));
        st->last_len = N;
    } else if (st->last_len != N) {
        uint8_t prev[513];
        float ratio;
        int i;

        memcpy(prev, st->last, sizeof(prev));
        ratio = st->last_len / (float)N;
        for (i = 0; i < N / 2; i++) {
            int k = (int)(i * ratio);
            st->last[i] = prev[k];
        }
        for (; i < (int)sizeof(st->last); i++)
            st->last[i] = 0;
        st->last_len = N;
    }

    return (uint32_t)st->last[0] << 12;
}

/* arith_get_context: compute the context of the current 2-tuple. */
static uint32_t aac_ac_get_context(aac_usac_ac_state* st, int i,
                                   uint16_t state_pre, uint8_t cur1,
                                   uint8_t cur2, uint8_t cur3,
                                   uint16_t* new_pre)
{
    uint32_t c = state_pre >> 8;
    c += (uint32_t)st->last[i + 1] << 8;
    c = (c << 4) + cur1;
    *new_pre = (uint16_t)c;

    if (i > 3 && (cur3 + cur2 + cur1) < 5)
        return c + 0x10000u;
    return c;
}

static void aac_ac_update_context(aac_usac_ac_state* st, int idx,
                                  int a, int b, uint8_t cur[4])
{
    int sum = a + b + 1;

    /* Cap in int precision: after long escape runs a + b + 1 can exceed 255,
     * so truncating to uint8_t first would wrap and corrupt the context. */
    if (sum > 0xF)
        sum = 0xF;
    cur[0] = (uint8_t)sum;
    cur[3] = cur[2];
    cur[2] = cur[1];
    cur[1] = cur[0];
    st->last[idx] = cur[0];
}

/* arith_finish_context: fill the remaining context with the stop value. */
static void aac_ac_finish_context(aac_usac_ac_state* st, int offset, int N)
{
    int i;

    for (i = offset; i < N / 2; i++)
        st->last[i] = 1;
    for (; i < (int)sizeof(st->last); i++)
        st->last[i] = 0;
}

int aac_usac_ac_decode(aac_usac_ac_state* st, aac_bits* b, float* coef,
                       int lg, int lg_max, int reset)
{
    uint16_t state_pre = 0;
    uint8_t cur[4];
    int i;
    size_t init_pos;

    if (lg_max < lg || lg_max > 2048)
        return -1;

    memset(coef, 0, sizeof(float) * (size_t)lg_max);

    aac_ac_map_process(st, reset, lg_max);
    cur[3] = 0;
    cur[2] = 0;
    cur[1] = 0;
    cur[0] = 1;
    state_pre = (uint16_t)((uint32_t)st->last[0] << 12);

    if (lg <= 0) {
        aac_ac_finish_context(st, 0, lg_max);
        return 0;
    }

    init_pos = b->bitpos;
    aac_ac_engine e;
    aac_ac_engine_init(&e, b);
#ifdef AAC_USAC_DEBUG
    fprintf(stderr, "ac: init pos=%zu val=%u reset=%d\n", init_pos,
            (unsigned)e.val, reset);
#endif

    for (i = 0; i < lg / 2; i++) {
        uint32_t c;
        int lvl, esc_nb, m;
        int a, bb;

        c = aac_ac_get_context(st, i, state_pre, cur[1], cur[2], cur[3],
                               &state_pre);

        for (lvl = esc_nb = 0;;) {
            int pki = aac_ac_get_pk(c + ((uint32_t)esc_nb << 17));
            m = aac_ac_decode(&e, b, aac_usac_ac_msb_cdfs[pki], 17);
            if (m < AAC_USAC_AC_ESCAPE)
                break;
            lvl++;
            if (lvl > 23)
                return -1;
            if (lvl > 7)
                esc_nb = 7;
            else
                esc_nb = lvl;
        }

        bb = m >> 2;
        a = m & 3;

        if (!m) {
            if (esc_nb)
                break;              /* stop symbol */
            a = 0;
            bb = 0;
        }

        /* least significant bits, interleaved a/b */
        for (int l = lvl; l > 0; l--) {
            int lsbidx = !a ? 1 : (!bb ? 0 : 2);
            int r = aac_ac_decode(&e, b, aac_usac_ac_lsb_cdfs[lsbidx], 4);
            a = (a << 1) | (r & 1);
            bb = (bb << 1) | ((r >> 1) & 1);
        }

        /* Dequantize: the USAC arithmetic spectrum stores |q|^(4/3) values,
         * the caller only applies the band scale factor. */
        coef[2 * i + 0] = (a >= 0) ? aac_pow43(a) : -aac_pow43(-a);
        coef[2 * i + 1] = (bb >= 0) ? aac_pow43(bb) : -aac_pow43(-bb);

        aac_ac_update_context(st, i, a, bb, cur);
    }

    /* The arithmetic decoder prefetches up to 14 bits past the actual data
     * during renormalization; rewind them so the signs are read from the
     * correct position. */
    if (b->bitpos >= init_pos + 14u)
        b->bitpos -= 14u;
    else
        b->bitpos = init_pos;
    if (b->bitpos > b->size * 8u) {
        b->bitpos = b->size * 8u;
        b->error = 1;
    }

    aac_ac_finish_context(st, i, lg_max);

#ifdef AAC_USAC_DEBUG
    fprintf(stderr, "ac: lg=%d lg_max=%d reset=%d stop_i=%d consumed=%zu\n",
            lg, lg_max, reset, i, b->bitpos - init_pos);
#endif

    for (; i < lg_max / 2; i++) {
        coef[2 * i + 0] = 0.0f;
        coef[2 * i + 1] = 0.0f;
    }

    for (i = 0; i < lg; i++) {
        if (coef[i] != 0.0f && !aac_bits_read(b, 1))
            coef[i] = -coef[i];
    }

    return b->error ? -1 : 0;
}
