/*
 * aacdec - Spectral Band Replication (HE-AAC) decoder.
 *
 * Clean-room implementation of the MPEG-4 SBR tools used by HE-AAC and
 * HE-AACv2: SBR header/grid/envelope/noise parsing, 32-band QMF analysis,
 * 64-band QMF synthesis, inverse filtering, high-frequency generation and
 * envelope/noise adjustment.
 */
#ifndef AAC_SBR_H
#define AAC_SBR_H

#include <stddef.h>
#include <stdint.h>

#include "aac_bits.h"
#include "aac_ps.h"

#define AAC_SBR_MAX_CHANNELS    2
#define AAC_SBR_MAX_BANDS       48
#define AAC_SBR_MAX_ENV         8
#define AAC_SBR_MAX_NOISE       3
#define AAC_SBR_MAX_PATCHES     6
#define AAC_SBR_MAX_SLOTS       40
#define AAC_SBR_MAX_HF_BANDS    64
#define AAC_SBR_MAX_TIME_SLOTS  32

/* The SBR element state is intentionally self contained so that it can be
 * embedded in the AAC core without dynamic allocation. */
typedef struct
{
    int amp_res;
    int start_freq;
    int stop_freq;
    int xover_band;
    int freq_scale;
    int alter_scale;
    int noise_bands;
    int limiter_bands;
    int limiter_gains;
    int interpol_freq;
    int smoothing_mode;
} aac_sbr_header;

/* Default SBR header parameters carried in UsacConfig (sbrConfig).  They
 * are applied when a USAC frame signals "default header". */
typedef struct
{
    int start_freq;
    int stop_freq;
    int freq_scale;
    int alter_scale;
    int noise_bands;
    int limiter_bands;
    int limiter_gains;
    int interpol_freq;
    int smoothing_mode;
} aac_sbr_usac_defaults;

typedef struct
{
    /* Bitstream data. */
    int frame_class;
    int num_env;
    int freq_res[9];            /* index 0 is previous frame's last value */
    int num_noise;
    int df_env[9];
    int df_noise[3];
    int invf_mode[2][5];        /* [new/old][band] */
    int add_harmonic_flag;
    int add_harmonic[48];
    int amp_res;
    int temp_shape[9];          /* USAC bs_interTes temp shape (parsed,
                                   not applied) */
    int temp_shape_mode[9];

    int t_env[9];
    int t_q[3];
    int t_env_num_env_old;
    int e_a[2];                 /* l_A and l_APrev */

    int env_facs_q[9][48];
    int noise_facs_q[3][5];
    float env_facs[9][48];
    float noise_facs[3][5];

    float bw_array[5];
    int f_indexnoise;
    int f_indexsine;

    /* Persistent QMF and synthesis state. */
    float ana_hist[320];
    float W[2][32][32][2];      /* [frame buffer][slot][band][re,im] */
    int Ypos;                   /* also selects the W buffer */
    float X_low[32][40][2];
    float X_high[64][40][2];
    float Y[2][38][64][2];      /* [frame buffer][slot][band][re,im] */
    float synth_buf[2048 + 640];

    /* Envelope/noise gain smoothing history. */
    float g_temp[42][48];
    float q_temp[42][48];
    int s_indexmapped[9][48];
} aac_sbr_data;

typedef struct
{
    int channels;               /* 1 for SCE, 2 for CPE */
    int started;                /* an SBR header has been seen */
    int active;                 /* this frame contains SBR data */
    int initialized;            /* frequency tables are valid */
    int ps_present;             /* PS extension data was found */
    int core_rate;
    int out_rate;
    int frame_size;
    int num_time_slots;

    int reset;

    aac_sbr_header hdr;
    int hdr_valid;

    int k[3];
    int kx[2];                  /* [previous, current] */
    int m[2];
    int n[2];
    int n_q;
    int n_master;
    int n_lim;

    uint16_t f_master[49];
    uint16_t f_tablelow[25];
    uint16_t f_tablehigh[49];
    uint16_t f_tablenoise[6];
    uint16_t f_tablelim[30];

    int num_patches;
    uint8_t patch_num_subbands[6];
    uint8_t patch_start_subband[6];

    int bs_coupling;
    int bs_amp_res_header;
    int bs_limiter_bands;
    int bs_limiter_gains;
    int bs_interpol_freq;
    int bs_smoothing_mode;
    int usac;                   /* USAC (AOT 42) SBR framing */
    int inter_tes;              /* USAC bs_interTes envelope temp shape */
    int xover_band_set;         /* USAC SbrInfo carried bs_xover_band */

    aac_sbr_data data[2];
    aac_ps_decoder ps;

    /* Shared per-frame scratch (channels are processed sequentially). */
    float alpha0[64][2];
    float alpha1[64][2];
    float e_origmapped[8][48];
    float q_mapped[8][48];
    float s_mapped[8][48];
    float e_curr[8][48];
    float q_m[8][48];
    float s_m[8][48];
    float gain[8][48];
} aac_sbr_decoder;

/* Initialize/reset an SBR decoder.  out_rate is normally 2*core_rate. */
int aac_sbr_open(aac_sbr_decoder* sbr, int channels, int core_rate,
                 int out_rate, int frame_size);
void aac_sbr_close(aac_sbr_decoder* sbr);
void aac_sbr_reset(aac_sbr_decoder* sbr);

/* Parse one SBR extension payload.  element_type is 0 for SCE and 1 for CPE.
 * len_bits is the number of bits available for the SBR data (excluding the
 * extension type).  crc indicates EXT_SBR_DATA_CRC. */
int aac_sbr_parse(aac_sbr_decoder* sbr, aac_bits* b, int len_bits,
                  int element_type, int crc);

/* Parse one USAC sbrData() payload (ISO/IEC 23003-3 7.3.5).  sbr_ch is the
 * number of SBR channels (1 = SCE or downmixed CPE, 2 = CPE).  indep_flag
 * reflects the UsacFrame independency flag.  dflt carries the sbrConfig
 * defaults from UsacConfig. */
int aac_sbr_parse_usac(aac_sbr_decoder* sbr, aac_bits* b, int sbr_ch,
                       int indep_flag, const aac_sbr_usac_defaults* dflt);

/* Apply SBR to one core frame.  core_pcm holds one pointer per channel and
 * out_pcm receives out_samples interleaved per channel plane.  When SBR data
 * is absent the decoder still performs the 2x QMF upsampling. */
int aac_sbr_apply(aac_sbr_decoder* sbr, const float* core_pcm[2],
                  int core_samples, float* out_pcm[2], int out_samples);

int aac_sbr_is_active(const aac_sbr_decoder* sbr);

#endif /* AAC_SBR_H */
