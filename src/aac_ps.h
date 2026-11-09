/*
 * aacdec - MPEG-4 Parametric Stereo (HE-AACv2).
 *
 * Baseline 20-band PS: IID/ICC parsing, hybrid analysis, all-pass
 * decorrelation and stereo mixing.  IPD/OPD and the 34-band hybrid mode are
 * mapped to the 20-band mode.
 */
#ifndef AAC_PS_H
#define AAC_PS_H

#include <stdint.h>

#include "aac_bits.h"

#define AAC_PS_MAX_ENV        5
#define AAC_PS_MAX_PAR        34
#define AAC_PS_MAX_SSB        91
#define AAC_PS_MAX_AP_BANDS   50
#define AAC_PS_MAX_DELAY      14
#define AAC_PS_AP_LINKS       3
#define AAC_PS_AP_DELAY       5
#define AAC_PS_QMF_SLOTS      32

typedef struct
{
    int start;
    int enable_iid;
    int iid_quant;
    int nr_iid_par;
    int iid_symbols;
    int freq_res_iid;
    int nr_ipdopd_par;
    int enable_icc;
    int icc_mode;
    int nr_icc_par;
    int icc_symbols;
    int freq_res_icc;
    int enable_ext;
    int enable_ipdopd;
    int frame_class;
    int num_env;
    int num_env_old;
    int border_position[AAC_PS_MAX_ENV + 1];
    int8_t iid_par[AAC_PS_MAX_ENV][AAC_PS_MAX_PAR];
    int8_t icc_par[AAC_PS_MAX_ENV][AAC_PS_MAX_PAR];
    int8_t ipd_par[AAC_PS_MAX_ENV][AAC_PS_MAX_PAR];
    int8_t opd_par[AAC_PS_MAX_ENV][AAC_PS_MAX_PAR];
    int is34bands;
    int is34bands_old;

    /* Hybrid analysis history. */
    float in_buf[5][44][2];

    /* Decorrelator history. */
    float delay[AAC_PS_MAX_SSB][AAC_PS_QMF_SLOTS + AAC_PS_MAX_DELAY][2];
    float ap_delay[AAC_PS_MAX_AP_BANDS][AAC_PS_AP_LINKS]
                  [AAC_PS_QMF_SLOTS + AAC_PS_AP_DELAY][2];
    float peak_decay_nrg[AAC_PS_MAX_PAR];
    float power_smooth[AAC_PS_MAX_PAR];
    float peak_decay_diff_smooth[AAC_PS_MAX_PAR];

    /* Mixing coefficients. */
    float H11[2][AAC_PS_MAX_ENV + 1][AAC_PS_MAX_PAR];
    float H12[2][AAC_PS_MAX_ENV + 1][AAC_PS_MAX_PAR];
    float H21[2][AAC_PS_MAX_ENV + 1][AAC_PS_MAX_PAR];
    float H22[2][AAC_PS_MAX_ENV + 1][AAC_PS_MAX_PAR];

    /* Hybrid subband buffers. */
    float Lbuf[AAC_PS_MAX_SSB][AAC_PS_QMF_SLOTS][2];
    float Rbuf[AAC_PS_MAX_SSB][AAC_PS_QMF_SLOTS][2];
    int8_t opd_hist[AAC_PS_MAX_PAR];
    int8_t ipd_hist[AAC_PS_MAX_PAR];
} aac_ps_decoder;

void aac_ps_reset(aac_ps_decoder* ps);

/* Parse one PS extension payload.  bits_left is the number of bits available
 * after the 2-bit extension id.  Returns the number of bits consumed. */
int aac_ps_parse(aac_ps_decoder* ps, aac_bits* b, int bits_left);

/* Apply PS to one SBR frame.  L is modified in place; R is produced.  top is
 * the highest active QMF band (kx + m). */
void aac_ps_apply(aac_ps_decoder* ps, float Lre[38][64],
                  float Lim[38][64], float Rre[38][64],
                  float Rim[38][64], int top);

#endif /* AAC_PS_H */
