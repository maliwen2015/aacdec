/*
 * aacdec - a self-contained MPEG-2 / MPEG-4 AAC-LC / HE-AAC / AAC-LD /
 * AAC-ELD decoder with a packet-based API.
 *
 * This is a clean-room implementation of the AAC-LC tools and the HE-AAC SBR
 * extension: ADTS, ADIF, LATM/LOAS and raw access units, mono to 7.1 channel
 * configurations (including PCE), M/S and intensity stereo, PNS, TNS, pulse
 * data, sine/KBD windows, the 1024/960-sample core filterbank and the
 * 32/64-band SBR QMF.  HE-AACv2 PS supports 20/34-band hybrid modes,
 * IID/ICC mixing and IPD/OPD phase parameters.  AAC-LD (AOT 23) and AAC-ELD
 * (AOT 39) add the 512/480-sample low-delay filterbank; LD-SBR, ELDv2, ER
 * protection (epConfig > 0) and MPEG-D USAC are not implemented by this core.
 * MPEG-4 dynamic range control is parsed and exported, and damaged access
 * units can be replaced by a faded copy of the previous frame.  ER AAC-LC
 * (AOT 17), ADTS header CRC checking and LATM AudioMuxVersion 0/1 with
 * multiple subframes are supported.
 *
 * Public API.  The interface deliberately mirrors the mpeg2dec / h264dec /
 * h265dec APIs: open a decoder, send packets, receive reference-counted
 * frames, flush, close.
 */
#ifndef AACDEC_H
#define AACDEC_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes */
#define AAC_DEC_ERROR       -1   /* decoding / transport error */
#define AAC_DEC_DEFAULT      0   /* no output available right now */
#define AAC_DEC_GOT_FRAME    1   /* got a decoded audio frame */
#define AAC_DEC_GOT_HEADER   2   /* got a new stream configuration */
#define AAC_DEC_END          3   /* end of stream after flush */

/* Transport formats.  AAC_TRANSPORT_AUTO probes ADTS, ADIF and LOAS/LATM from
 * the first bytes; when out-of-band config bytes are supplied it defaults to
 * raw access units. */
enum {
    AAC_TRANSPORT_AUTO = 0,
    AAC_TRANSPORT_RAW,          /* raw access units, config (ASC) required */
    AAC_TRANSPORT_ADTS,         /* ADTS sync stream (.aac) */
    AAC_TRANSPORT_ADIF,         /* ADIF stream */
    AAC_TRANSPORT_LATM,         /* LATM/LOAS sync stream (in-band config) */
    AAC_TRANSPORT_LOAS = AAC_TRANSPORT_LATM,
    AAC_TRANSPORT_LATM_MCP0,    /* packetized AudioMuxElement, out-of-band SMC */
    AAC_TRANSPORT_LATM_MCP1     /* packetized AudioMuxElement, in-band config */
};

/* PCM output formats */
enum {
    AAC_OUTPUT_S16 = 0,         /* interleaved signed 16-bit (native) */
    AAC_OUTPUT_S16_PLANAR,      /* one S16 plane per channel */
    AAC_OUTPUT_F32,             /* interleaved 32-bit float */
    AAC_OUTPUT_F32_PLANAR       /* one float plane per channel */
};

/* Decoder parameter flags */
#define AAC_DEC_FLAG_DISABLE_CONCEALMENT 0x1

/* Channel ordering of the output buffer */
enum {
    AAC_CHANNEL_ORDER_WAV = 0,  /* WAV/MPEG-4 order (default) */
    AAC_CHANNEL_ORDER_MPEG      /* PCE / ISO 13818-7 canonical order */
};

/* Audio object types (ISO/IEC 14496-3) */
#define AAC_AOT_AAC_MAIN        1
#define AAC_AOT_AAC_LC          2
#define AAC_AOT_AAC_SSR         3
#define AAC_AOT_AAC_LTP         4
#define AAC_AOT_SBR             5
#define AAC_AOT_ER_AAC_LC       17
#define AAC_AOT_ER_AAC_LD       23
#define AAC_AOT_PS              29
#define AAC_AOT_ER_AAC_ELD      39
#define AAC_AOT_USAC            42

/* Speaker position tags, compatible with the MPEG-4 channel mapping. */
#define AAC_CHANNEL_TYPE_NONE          0x00
#define AAC_CHANNEL_TYPE_FRONT         0x01
#define AAC_CHANNEL_TYPE_SIDE          0x02
#define AAC_CHANNEL_TYPE_BACK          0x03
#define AAC_CHANNEL_TYPE_LFE           0x04
#define AAC_CHANNEL_TYPE_TOP           0x10
#define AAC_CHANNEL_TYPE_FRONT_TOP     0x11
#define AAC_CHANNEL_TYPE_SIDE_TOP      0x12
#define AAC_CHANNEL_TYPE_BACK_TOP      0x13
#define AAC_CHANNEL_TYPE_BOTTOM        0x20
#define AAC_CHANNEL_TYPE_FRONT_BOTTOM  0x21
#define AAC_CHANNEL_TYPE_SIDE_BOTTOM   0x22
#define AAC_CHANNEL_TYPE_BACK_BOTTOM   0x23

#define AAC_DEC_MAX_CHANNELS   8
#define AAC_DEC_MAX_FRAME_SIZE 4096

/* Stream information, updated whenever the configuration changes. */
typedef struct aac_stream_info
{
    int sample_rate;             /* output sample rate in Hz */
    int frame_size;              /* PCM samples per channel in one frame */
    int channels;                /* number of output channels */

    int aot;                     /* audio object type (AAC_AOT_*) */
    int ext_aot;                 /* extension object type, 0 if none */
    int ext_sample_rate;         /* SBR/PS output rate, 0 if none */
    int channel_config;          /* ISO channelConfiguration, 0 = PCE */
    int profile;                 /* MPEG-2 profile, -1 for MPEG-4 */
    int bitrate;                 /* instantaneous bitrate in bit/s */

    int aac_sample_rate;         /* AAC core sample rate */
    int aac_frame_size;          /* AAC core frame length */
    int aac_channels;            /* channels before PS/MPS */
    int ep_config;               /* error protection config, -1 = none */
    unsigned output_delay;       /* additional decoder delay in samples */

    int channel_types[AAC_DEC_MAX_CHANNELS];
    int channel_indices[AAC_DEC_MAX_CHANNELS];

    int num_lost_access_units;   /* estimated lost AUs after a sync error */
    int64_t num_total_bytes;
    int64_t num_bad_bytes;
    int64_t num_total_access_units;
    int64_t num_bad_access_units;

    int drc_program_ref_level;   /* -1 when absent */
    int drc_presentation_mode;   /* -1 when absent */
    int output_loudness;         /* -1 when absent */
    unsigned flags;              /* internal codec flags (AC_*) */
} aac_stream_info;

/* Input packet.  The decoder buffers data internally, so a packet may contain
 * any number of bytes: whole access units, partial access units, or several
 * access units.  Raw access-unit transports (AAC_TRANSPORT_RAW and the
 * packetized LATM modes) require one complete access unit per packet. */
typedef struct aac_packet_t
{
    const uint8_t* data;
    int len;
    int64_t pts;
    int64_t dts;
    void* opaque;
} aac_packet_t;

/* Decoded audio frame.  The decoder owns the sample buffer until
 * aac_dec_frame_unref() is called.
 *
 * For interleaved output num_planes is 1 and planes[0] holds
 * channels interleaved samples; strides[0] is the byte distance between two
 * consecutive sample frames (channels * bytes_per_sample).
 * For planar output num_planes equals channels, planes[c] points at channel c
 * and strides[c] is bytes_per_sample. */
typedef struct aac_frame_t
{
    uint8_t* planes[AAC_DEC_MAX_CHANNELS];
    int strides[AAC_DEC_MAX_CHANNELS];
    int num_planes;

    int channels;
    int frame_size;             /* samples per channel */
    int sample_rate;
    int format;                 /* AAC_OUTPUT_* */
    int bytes_per_sample;       /* 2 or 4 */
    int is_planar;

    int error;                  /* 0 = clean, 1 = concealed/corrupt frame */

    int64_t pts;
    int64_t dts;
    void* opaque;

    void* dec_frame;            /* internal handle, pass to frame_unref */
} aac_frame_t;

/* Decoder parameters.  A zero-initialised structure selects the defaults:
 * automatic transport detection, interleaved S16 output and WAV channel
 * order. */
typedef struct aac_param_t
{
    int transport;               /* AAC_TRANSPORT_* */
    int output_format;           /* AAC_OUTPUT_* */
    const uint8_t* config;       /* AudioSpecificConfig / StreamMuxConfig */
    int config_len;              /* bytes in config, 0 = none */
    int layers;                  /* reserved, currently ignored */
    int channel_order;           /* AAC_CHANNEL_ORDER_* */
    int flags;                   /* AAC_DEC_FLAG_* */
} aac_param_t;

#ifdef AACDECLIB_SHARED
#if defined(_WIN32)
#if defined(AACDECLIB_BUILD)
#define AACDEC_EXPORT __declspec(dllexport)
#else
#define AACDEC_EXPORT __declspec(dllimport)
#endif
#else
#define AACDEC_EXPORT __attribute__((visibility("default")))
#endif
#else
#define AACDEC_EXPORT
#endif

/* Open a decoder.  Returns NULL on invalid parameters or out of memory. */
AACDEC_EXPORT void* aac_dec_open(const aac_param_t* param);

/* Supply or replace the out-of-band configuration (ASC/SMC).  This must be
 * done before the first access unit is decoded.  Returns 0 on success. */
AACDEC_EXPORT int aac_dec_set_config(void* decoder, const void* config, int len);

/* Queue compressed bytes.  Returns 0 on success.  The packet data is copied;
 * the caller may free it as soon as the call returns. */
AACDEC_EXPORT int aac_dec_send_packet(void* decoder, const aac_packet_t* packet);

/* Fetch a configuration update or a decoded frame.  On AAC_DEC_GOT_HEADER a
 * frame decoded from the same access unit is retained internally and returned
 * by the next call. */
AACDEC_EXPORT int aac_dec_recv_frame(void* decoder, aac_stream_info* info,
                                     aac_frame_t* frame);

/* Drain all buffered access units at the end of the stream.  Returns
 * AAC_DEC_END when no more frames are available. */
AACDEC_EXPORT int aac_dec_flush(void* decoder, aac_stream_info* info,
                                aac_frame_t* frame);

/* Generate one concealment frame after the caller detects a lost access
 * unit.  The frame repeats and fades the last successfully decoded frame and
 * is returned with frame->error = 1.  Call this between recv calls, after the
 * decoder has drained all buffered access units.  Returns AAC_DEC_GOT_FRAME
 * or a negative error when no concealment history is available. */
AACDEC_EXPORT int aac_dec_conceal_frame(void* decoder, aac_stream_info* info,
                                        aac_frame_t* frame);

/* Release a frame returned by recv_frame/flush. */
AACDEC_EXPORT void aac_dec_frame_unref(void* decoder, aac_frame_t* frame);

/* Copy the most recent stream information.  Returns 0 on success, -1 when no
 * configuration has been seen yet. */
AACDEC_EXPORT int aac_dec_get_stream_info(void* decoder, aac_stream_info* info);

/* Destroy the decoder and free all buffers. */
AACDEC_EXPORT void aac_dec_close(void* decoder);

/* Version and error strings. */
AACDEC_EXPORT const char* aac_dec_version(void);
AACDEC_EXPORT const char* aac_dec_strerror(int code);

#ifdef __cplusplus
}
#endif

#endif /* AACDEC_H */
