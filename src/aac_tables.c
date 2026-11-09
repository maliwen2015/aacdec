/*
 * aacdec - standard AAC tables and window generation.
 */
#include "aac_tables.h"
#include "aac_ld_tables.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "aac_tables_data.h"

const int aac_sample_rates[AAC_SR_INDEX_MAX] = {
    96000, 88200, 64000, 48000, 44100, 32000,
    24000, 22050, 16000, 12000, 11025, 8000, 7350
};

static const uint16_t* const swb_long_1024[AAC_SR_INDEX_MAX] = {
    aac_swb_1024_96, aac_swb_1024_96, aac_swb_1024_64,
    aac_swb_1024_48, aac_swb_1024_48, aac_swb_1024_32,
    aac_swb_1024_24, aac_swb_1024_24, aac_swb_1024_16,
    aac_swb_1024_16, aac_swb_1024_16, aac_swb_1024_8,
    aac_swb_1024_8
};

static const uint16_t* const swb_short_128[AAC_SR_INDEX_MAX] = {
    aac_swb_128_96, aac_swb_128_96, aac_swb_128_96,
    aac_swb_128_48, aac_swb_128_48, aac_swb_128_48,
    aac_swb_128_24, aac_swb_128_24, aac_swb_128_16,
    aac_swb_128_16, aac_swb_128_16, aac_swb_128_8,
    aac_swb_128_8
};

static const uint16_t* const swb_long_960[AAC_SR_INDEX_MAX] = {
    aac_swb_960_96, aac_swb_960_96, aac_swb_960_64,
    aac_swb_960_48, aac_swb_960_48, aac_swb_960_48,
    aac_swb_960_24, aac_swb_960_24, aac_swb_960_16,
    aac_swb_960_16, aac_swb_960_16, aac_swb_960_8,
    aac_swb_960_8
};

static const uint16_t* const swb_short_120[AAC_SR_INDEX_MAX] = {
    aac_swb_120_96, aac_swb_120_96, aac_swb_120_96,
    aac_swb_120_48, aac_swb_120_48, aac_swb_120_48,
    aac_swb_120_24, aac_swb_120_24, aac_swb_120_16,
    aac_swb_120_16, aac_swb_120_16, aac_swb_120_8,
    aac_swb_120_8
};

typedef struct
{
    float long_half[1024];
    float long_kbd[1024];
    float short_half[128];
    float short_kbd[128];
    int frame_size;
} aac_window_set;

static aac_window_set win_1024;
static aac_window_set win_960;
static aac_window_set win_512;
static aac_window_set win_480;
static float pow2_quarter[512];
static float pow43[8192];
static int tables_ready;

static double bessel_i0(double x)
{
    double sum = 1.0;
    double term = 1.0;
    double x2 = x * x * 0.25;
    int k;

    for (k = 1; k < 128; k++) {
        term *= x2 / ((double)k * (double)k);
        sum += term;
        if (term < 1e-18 * sum)
            break;
    }
    return sum;
}

static void kbd_half(float* out, int half, double alpha)
{
    const int n = half * 2;
    int i;
    double* temp;
    double total = 0.0;
    double sum = 0.0;

    temp = (double*)malloc(((size_t)half + 1) * sizeof(*temp));
    if (temp == NULL)
        return;

    for (i = 0; i <= half; i++) {
        double x = 4.0 * (double)i / (double)n - 1.0;
        double arg = 3.14159265358979323846 * alpha *
                     sqrt(1.0 - x * x);
        temp[i] = bessel_i0(arg);
        total += temp[i];
    }

    for (i = 0; i < half; i++) {
        sum += temp[i];
        out[i] = (float)sqrt(sum / total);
    }
    free(temp);
}

static void window_set_init(aac_window_set* w, int frame_size)
{
    int i;
    int short_len = frame_size / 8;

    w->frame_size = frame_size;

    for (i = 0; i < frame_size; i++)
        w->long_half[i] =
            (float)sin(3.14159265358979323846 / (2.0 * frame_size) *
                       ((double)i + 0.5));
    for (i = 0; i < short_len; i++)
        w->short_half[i] =
            (float)sin(3.14159265358979323846 / (2.0 * short_len) *
                       ((double)i + 0.5));

    kbd_half(w->long_kbd, frame_size, 4.0);
    kbd_half(w->short_kbd, short_len, 6.0);
}

static const uint16_t* const swb_512[AAC_SR_INDEX_MAX] = {
    NULL, NULL, NULL,
    aac_ld_swb512_48, aac_ld_swb512_48, aac_ld_swb512_32,
    aac_ld_swb512_24, aac_ld_swb512_24, NULL,
    NULL, NULL, NULL, NULL
};

static const uint16_t* const swb_480[AAC_SR_INDEX_MAX] = {
    NULL, NULL, NULL,
    aac_ld_swb480_48, aac_ld_swb480_48, aac_ld_swb480_32,
    aac_ld_swb480_24, aac_ld_swb480_24, NULL,
    NULL, NULL, NULL, NULL
};

const uint16_t* aac_swb_offsets(int frame_size, int short_window, int sr_index)
{
    if (sr_index < 0 || sr_index >= AAC_SR_INDEX_MAX)
        sr_index = 4;

    if (frame_size == 512)
        return short_window ? NULL : swb_512[sr_index];
    if (frame_size == 480)
        return short_window ? NULL : swb_480[sr_index];
    if (frame_size == 960)
        return short_window ? swb_short_120[sr_index] : swb_long_960[sr_index];
    return short_window ? swb_short_128[sr_index] : swb_long_1024[sr_index];
}

int aac_num_swb(int frame_size, int short_window, int sr_index)
{
    if (sr_index < 0 || sr_index >= AAC_SR_INDEX_MAX)
        sr_index = 4;
    if (frame_size == 512)
        return short_window ? 0 : aac_ld_num_swb_512[sr_index];
    if (frame_size == 480)
        return short_window ? 0 : aac_ld_num_swb_480[sr_index];
    if (frame_size == 960)
        return short_window ? aac_num_swb_120[sr_index]
                            : aac_num_swb_960[sr_index];
    return short_window ? aac_num_swb_128[sr_index]
                        : aac_num_swb_1024[sr_index];
}

int aac_tns_max_bands(int frame_size, int short_window, int sr_index)
{
    if (sr_index < 0 || sr_index >= AAC_SR_INDEX_MAX)
        sr_index = 4;
    if (frame_size == 512)
        return short_window ? 0 : aac_ld_tns_512[sr_index];
    if (frame_size == 480)
        return short_window ? 0 : aac_ld_tns_480[sr_index];
    return short_window ? aac_tns_max_bands_128[sr_index]
                        : aac_tns_max_bands_1024[sr_index];
}

const float* aac_window_half(int frame_size, int short_window, int shape)
{
    aac_window_set* w;
    if (frame_size == 512)
        w = &win_512;
    else if (frame_size == 480)
        w = &win_480;
    else if (frame_size == 960)
        w = &win_960;
    else
        w = &win_1024;
    if (shape)
        return short_window ? w->short_kbd : w->long_kbd;
    return short_window ? w->short_half : w->long_half;
}

float aac_pow2_quarter(int e)
{
    if (e < -256)
        e = -256;
    if (e > 255)
        e = 255;
    return pow2_quarter[e + 256];
}

float aac_pow43(int q)
{
    if (q < 0)
        q = 0;
    if (q > 8191)
        q = 8191;
    return pow43[q];
}

int aac_tables_init(void)
{
    int i;

    if (tables_ready)
        return 0;

    window_set_init(&win_1024, 1024);
    window_set_init(&win_960, 960);
    window_set_init(&win_512, 512);
    window_set_init(&win_480, 480);

    for (i = 0; i < 512; i++)
        pow2_quarter[i] = (float)exp2(((double)(i - 256)) * 0.25);

    pow43[0] = 0.0f;
    for (i = 1; i < 8192; i++)
        pow43[i] = (float)pow((double)i, 4.0 / 3.0);

    tables_ready = 1;
    return 0;
}
