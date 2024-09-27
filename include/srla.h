#ifndef SRLA_H_INCLUDED
#define SRLA_H_INCLUDED

#include "srla_stdint.h"

/* Format version */
#define SRLA_FORMAT_VERSION         5

/* Codec version */
#define SRLA_CODEC_VERSION          8

/* Header size */
#define SRLA_HEADER_SIZE            29

/* Maximum number of channels that can be processed */
#define SRLA_MAX_NUM_CHANNELS       8

/* Maximum coefficient size */
#define SRLA_MAX_COEFFICIENT_ORDER  256

/* Number of parameter presets */
#define SRLA_NUM_PARAMETER_PRESETS  5


/* API result type */
typedef enum SRLAApiResultTag {
    SRLA_APIRESULT_OK = 0,                  /* success */
    SRLA_APIRESULT_INVALID_ARGUMENT,        /* Invalid argument */
    SRLA_APIRESULT_INVALID_FORMAT,          /* Invalid format */
    SRLA_APIRESULT_INSUFFICIENT_BUFFER,     /* Buffer size is insufficient */
    SRLA_APIRESULT_INSUFFICIENT_DATA,       /* Not enough data */
    SRLA_APIRESULT_PARAMETER_NOT_SET,       /* No parameters set */
    SRLA_APIRESULT_DETECT_DATA_CORRUPTION,  /* Data corruption detected */
    SRLA_APIRESULT_NG                       /* Unclassifiable failure */
} SRLAApiResult;

/* Header information */
struct SRLAHeader {
    uint32_t format_version;                        /* Format version */
    uint32_t codec_version;                         /* Encoder version */
    uint16_t num_channels;                          /* Number of channels */
    uint32_t num_samples;                           /* Total number of samples per channel */
    uint32_t sampling_rate;                         /* Sampling rate */
    uint16_t bits_per_sample;                       /* Number of bits per sample */
    uint32_t max_num_samples_per_block;             /* Maximum number of samples per block */
    uint8_t preset;                                 /* Parameter presets */
};

#endif /* SRLA_H_INCLUDED */
