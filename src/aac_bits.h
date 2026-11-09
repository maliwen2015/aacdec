/*
 * aacdec - bit reader used by the AAC syntax parser.
 */
#ifndef AAC_BITS_H
#define AAC_BITS_H

#include <stddef.h>
#include <stdint.h>

typedef struct aac_bits
{
    const uint8_t* data;
    size_t size;               /* size in bytes */
    size_t bitpos;             /* current bit position */
    int error;                 /* set when reading past the end */
} aac_bits;

static inline void aac_bits_init(aac_bits* b, const uint8_t* data, size_t size)
{
    b->data = data;
    b->size = size;
    b->bitpos = 0;
    b->error = 0;
}

static inline size_t aac_bits_left(const aac_bits* b)
{
    size_t total = b->size * 8u;
    return (b->bitpos < total) ? (total - b->bitpos) : 0u;
}

static inline uint32_t aac_bits_read(aac_bits* b, int n)
{
    uint32_t v = 0;

    if (n <= 0)
        return 0;
    if (n > 32)
        n = 32;

    while (n > 0) {
        size_t byte;
        int off, avail, take;
        uint8_t value;
        uint32_t mask;

        if (b->bitpos >= b->size * 8u) {
            b->error = 1;
            return 0;
        }
        byte = b->bitpos >> 3;
        off = (int)(b->bitpos & 7u);
        avail = 8 - off;
        take = (n < avail) ? n : avail;
        value = b->data[byte];
        mask = (1u << take) - 1u;
        v = (v << take) | ((uint32_t)(value >> (avail - take)) & mask);
        b->bitpos += (size_t)take;
        n -= take;
    }
    return v;
}

static inline uint32_t aac_bits_peek(aac_bits* b, int n)
{
    size_t save = b->bitpos;
    int save_error = b->error;
    uint32_t v = aac_bits_read(b, n);
    b->bitpos = save;
    b->error = save_error;
    return v;
}

static inline void aac_bits_skip(aac_bits* b, size_t n)
{
    b->bitpos += n;
    if (b->bitpos > b->size * 8u) {
        b->bitpos = b->size * 8u;
        b->error = 1;
    }
}

static inline void aac_bits_align(aac_bits* b)
{
    b->bitpos = (b->bitpos + 7u) & ~(size_t)7u;
    if (b->bitpos > b->size * 8u) {
        b->bitpos = b->size * 8u;
        b->error = 1;
    }
}

static inline size_t aac_bits_bytes_consumed(const aac_bits* b)
{
    return (b->bitpos + 7u) >> 3;
}

#endif /* AAC_BITS_H */
