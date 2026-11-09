/*
 * aacdec - USAC (xHE-AAC, ISO/IEC 23003-3) configuration parser.
 */
#ifndef AAC_USAC_H
#define AAC_USAC_H

#include <stddef.h>
#include <stdint.h>

#include "aac_bits.h"

#define AAC_USAC_MAX_ELEMENTS 32

#define AAC_USAC_ID_SCE 0
#define AAC_USAC_ID_CPE 1
#define AAC_USAC_ID_LFE 2
#define AAC_USAC_ID_EXT 3

/* usacExtElementType */
#define AAC_USAC_EXT_FILL          0
#define AAC_USAC_EXT_MPEGS         1
#define AAC_USAC_EXT_SAOC          2
#define AAC_USAC_EXT_ELEVATION     3
#define AAC_USAC_EXT_HOLEFILLER    4
#define AAC_USAC_EXT_AUDIOPREROLL  5
#define AAC_USAC_EXT_UNI_DRC       6

typedef struct
{
    int type;                   /* AAC_USAC_ID_* */

    /* UsacCoreConfig (SCE/CPE/LFE) */
    int tw_mdct;
    int noise_fill;
    int sbr_ratio;              /* 0 none, 1 = 4/3, 2 = 8/3, 3 = 2x */

    /* UsacChannelPairElementConf */
    int stereo_config_index;

    /* UsacExtElementConfig */
    int ext_type;
    int ext_default_len;
    int ext_payload_frag;

    /* sbrConfig() defaults (dflt_* values used to build the USAC SBR
     * header when the frame signals the default header). */
    int sbr_start_freq;
    int sbr_stop_freq;
    int sbr_freq_scale;
    int sbr_alter_scale;
    int sbr_noise_bands;
    int sbr_limiter_bands;
    int sbr_limiter_gains;
    int sbr_interpol_freq;
    int sbr_smoothing_mode;
} aac_usac_elem_config;

typedef struct
{
    int sample_rate;            /* usacSamplingFrequency */
    int core_sbr_frame_len_idx; /* 0..5 */
    int core_frame_len;         /* 768 or 1024 */
    int sbr_ratio;              /* 0 none, 1 = 4/3, 2 = 8/3, 3 = 2x */
    int channels;               /* from channelConfigurationIndex / config */
    int channel_config_index;   /* 0 = UsacChannelConfig() present */
    int nb_elements;            /* number of elements in UsacDecoderConfig */
    int sbr_present;            /* any element carries SBR */

    aac_usac_elem_config elems[AAC_USAC_MAX_ELEMENTS];
} aac_usac_config;

typedef struct
{
    aac_usac_config config;
    size_t config_bits;         /* bits consumed by UsacConfig() */
    int asc_end_bit;            /* bit position after the config (unaligned;
                                   LATM carries frameLengthType right after) */
} aac_usac_asc;

/* Parse AudioSpecificConfig() with audioObjectType == 42 (USAC).  Returns 0
 * on success and fills *asc.  The bit position after the config is stored in
 * asc->asc_end_bit so the caller can byte-align the remaining payload. */
int aac_usac_parse_asc(const uint8_t* data, size_t len, aac_usac_asc* asc);

/* Parse UsacConfig() starting at the current bit position of *b (which must
 * already be positioned after the audioObjectType escape bits).  On success
 * returns 0, *asc holds the config; asc->config_bits counts the bits
 * consumed by UsacConfig() itself and asc->asc_end_bit the absolute bit
 * position after byte alignment.  *b is advanced past the config. */
int aac_usac_parse_asc_bits_state(aac_bits* b, aac_usac_asc* asc);

#endif /* AAC_USAC_H */
