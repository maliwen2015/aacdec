/*
 * aacdec - standard AAC tables and window generation.
 */
#ifndef AAC_TABLES_H
#define AAC_TABLES_H

#include <stdint.h>

#define AAC_SR_INDEX_MAX 13

extern const int aac_sample_rates[AAC_SR_INDEX_MAX];

/* Returns a pointer to the scalefactor band offsets.  The returned array has
 * aac_num_swb() + 1 entries. */
const uint16_t* aac_swb_offsets(int frame_size, int short_window, int sr_index);
int aac_num_swb(int frame_size, int short_window, int sr_index);
int aac_tns_max_bands(int frame_size, int short_window, int sr_index);

/* Half-window tables.  frame_size is 1024 or 960; short_window selects the
 * 128/120-sample short window; shape is 0 = sine, 1 = KBD.  The returned
 * pointer refers to frame_size (or frame_size/8) floats. */
const float* aac_window_half(int frame_size, int short_window, int shape);

/* 2^(e/4), e in [-256, 255]. */
float aac_pow2_quarter(int e);

/* |q|^(4/3), q in [0, 8191]. */
float aac_pow43(int q);

int aac_tables_init(void);

#endif /* AAC_TABLES_H */
