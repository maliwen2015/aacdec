/*
 * aacdec - clean-room AAC-LC / AAC-LD / AAC-ELD decoder core.
 */
#ifndef AAC_CORE_H
#define AAC_CORE_H

#include <stddef.h>
#include <stdint.h>

#include "aac_sbr.h"
#include "aac_usac.h"
#include "aac_usac_ac.h"

#define AAC_CORE_MAX_CHANNELS 8
#define AAC_CORE_MAX_FRAME_SIZE 2048
#define AAC_CORE_MAX_SFB      64
#define AAC_CORE_MAX_GROUPS   8
#define AAC_CORE_MAX_WINDOWS  8
#define AAC_CORE_MAX_TNS_ORDER 20
#define AAC_CORE_MAX_TNS_FILTERS 3
#define AAC_CORE_MAX_PULSES   4

#define AAC_CORE_OK        0
#define AAC_CORE_NEED_MORE (-1)
#define AAC_CORE_ERROR     (-2)

typedef struct
{
    int aot;                    /* audio object type from ASC/ADTS */
    int core_aot;               /* core object type (2 for HE-AAC) */
    int profile;                /* MPEG-2 profile, -1 for MPEG-4 */
    int mpeg2;

    int sample_rate;            /* core sampling rate */
    int sample_rate_index;
    int ext_sample_rate;        /* SBR/PS output rate, 0 if none */
    int ext_aot;
    int has_sbr;
    int has_ps;

    int channel_config;
    int frame_length_flag;      /* 1 = 960 samples, 0 = 1024 */
    int frame_size;
    int short_len;              /* frame_size / 8 */
    int ep_config;

    int channels;
    int channel_types[AAC_CORE_MAX_CHANNELS];
    int channel_indices[AAC_CORE_MAX_CHANNELS];
    int output_order[AAC_CORE_MAX_CHANNELS]; /* WAV order -> internal channel */

    int pce_valid;
    int pce_element_count;
    struct
    {
        int type;               /* 0 = SCE, 1 = CPE, 3 = LFE */
        int is_cpe;
        int tag;
        int channel_type[2];
        int channel_index[2];
    } pce_elements[AAC_CORE_MAX_CHANNELS + 2];

    /* MPEG-4 dynamic range control (FIL extension payload type 11). */
    int drc_present;
    int drc_pce_instance_tag;
    int drc_prog_ref_present;
    int drc_prog_ref_level;
    int drc_band_incr;
    int drc_interpolation_scheme;
    int drc_num_bands;
    int drc_dyn_rng_sgn[17];
    int drc_dyn_rng_ctl[17];
} aac_config;

typedef struct
{
    float overlap[1024];
    int prev_window_shape;
    int initialized;
    float eld_z[3][2048];       /* AAC-ELD overlap-add history */
} aac_channel_state;

typedef struct
{
    int window_sequence;
    int window_shape;
    int max_sfb;
    int num_windows;
    int num_groups;
    int group_len[AAC_CORE_MAX_GROUPS];
    int num_swb;
    int tns_max_bands;
    const uint16_t* swb_offset;

    int global_gain;
    int band_type[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int sf[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int noise_energy[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int is_position[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int ms_used[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int ms_mask_present;

    int tns_present;
    int tns_nfilt[AAC_CORE_MAX_WINDOWS];
    int tns_length[AAC_CORE_MAX_WINDOWS][AAC_CORE_MAX_TNS_FILTERS];
    int tns_order[AAC_CORE_MAX_WINDOWS][AAC_CORE_MAX_TNS_FILTERS];
    int tns_direction[AAC_CORE_MAX_WINDOWS][AAC_CORE_MAX_TNS_FILTERS];
    float tns_coef[AAC_CORE_MAX_WINDOWS][AAC_CORE_MAX_TNS_FILTERS]
                  [AAC_CORE_MAX_TNS_ORDER];

    int pulse_present;
    int num_pulse;
    int pulse_pos[AAC_CORE_MAX_PULSES];
    int pulse_amp[AAC_CORE_MAX_PULSES];
} aac_ics;

typedef struct
{
    aac_ics ics;
    int q[1024];
    float spec[1024];
    float pcm[1024];
    int noise_used[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int intensity_used[AAC_CORE_MAX_GROUPS][AAC_CORE_MAX_SFB];
    int element_type;
    int element_tag;
    int channel_type;
    int channel_index;
    int usac_noise_level;       /* USAC noiseFillingData level (0 = off) */
    int usac_noise_offset;      /* USAC noiseFillingData offset */
} aac_channel_work;

typedef struct
{
    int channels;
    int channel_types[AAC_CORE_MAX_CHANNELS];
    int channel_indices[AAC_CORE_MAX_CHANNELS];
    int output_order[AAC_CORE_MAX_CHANNELS];
    int frame_size;
    int sample_rate;

    int element_count;
    int element_type[16];
    int element_tag[16];
    int element_first_channel[16];
    int element_channels[16];
    int element_sbr[16];        /* SBR payload was present for this element */

    aac_channel_work ch[AAC_CORE_MAX_CHANNELS];
} aac_frame_work;

typedef struct
{
    aac_config cfg;
    aac_channel_state state[AAC_CORE_MAX_CHANNELS];
    uint32_t pns_seed;
    int configured;
    int frames;
    int mpeg_order;             /* 1 = output canonical MPEG order */
    aac_frame_work work;

    /* USAC (AOT 42) decoder state. */
    int usac_valid;
    aac_usac_asc usac;
    uint32_t usac_noise_seed[AAC_CORE_MAX_CHANNELS];
    aac_usac_ac_state usac_ac[AAC_CORE_MAX_CHANNELS];

    /* Complex stereo prediction persistence (ISO/IEC 23003-3 7.7.2.3.2). */
    int usac_cplx_valid;
    int usac_cplx_win_groups_prev;
    int usac_cplx_win_seq_prev;
    int usac_cplx_win_shape_prev;
    int usac_cplx_win_groups_cur;
    int usac_cplx_win_seq_cur;
    int usac_cplx_win_shape_cur;
    int usac_cplx_alpha_re_prev[AAC_CORE_MAX_GROUPS * AAC_CORE_MAX_SFB];
    int usac_cplx_alpha_im_prev[AAC_CORE_MAX_GROUPS * AAC_CORE_MAX_SFB];
    float usac_cplx_prev_l[AAC_CORE_MAX_FRAME_SIZE];
    float usac_cplx_prev_r[AAC_CORE_MAX_FRAME_SIZE];
    int usac_cplx_prev_len;

    /* SBR state is kept per element and indexed by the element's first
     * output channel.  It is allocated on demand. */
    aac_sbr_decoder* sbr[AAC_CORE_MAX_CHANNELS];
    float pcm_hi[AAC_CORE_MAX_CHANNELS][AAC_CORE_MAX_FRAME_SIZE];
} aac_core;

int aac_core_init(aac_core* core);
void aac_core_close(aac_core* core);
void aac_core_reset(aac_core* core);

int aac_core_set_asc(aac_core* core, const uint8_t* asc, size_t len);
int aac_core_parse_asc_bits(aac_core* core, void* bits_ctx);
int aac_core_set_adts(aac_core* core, int profile, int sf_index,
                      int channel_config, int mpeg2);
int aac_core_set_adif_pce(aac_core* core, int profile, int sf_index,
                          const uint8_t* pce, size_t pce_len);
int aac_core_parse_pce_bits(aac_core* core, void* bits_ctx);
int aac_core_parse_dynamic_range(aac_core* core, void* bits_ctx);
int aac_core_apply_config(aac_core* core, int reset);

/* Decode one raw_data_block.  On success returns AAC_CORE_OK and sets
 * *consumed to the byte length, *channels to the PCM channel count and writes
 * interleaved float PCM in WAV channel order into pcm (frame_size * channels
 * samples). */
int aac_core_decode(aac_core* core, const uint8_t* data, size_t len,
                    size_t* consumed, float* pcm, int* channels,
                    int* sample_rate);

int aac_core_frame_size(const aac_core* core);

/* Frame size and sample rate of the decoder output.  These differ from the
 * AAC core values when SBR (HE-AAC) is active. */
int aac_core_output_frame_size(const aac_core* core);
int aac_core_output_sample_rate(const aac_core* core);

#endif /* AAC_CORE_H */
