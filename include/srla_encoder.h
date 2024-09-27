#ifndef SRLA_ENCODER_H_INCLUDED
#define SRLA_ENCODER_H_INCLUDED

#include "srla.h"
#include "srla_stdint.h"

/* Encoding parameters */
struct SRLAEncodeParameter {
    uint16_t num_channels; /* Number of channels of the input waveform */
    uint16_t bits_per_sample; /* Number of bits per sample of the input waveform */
    uint32_t sampling_rate; /* Sampling rate of input waveform */
    uint32_t min_num_samples_per_block; /* Minimum number of samples per block */
    uint32_t max_num_samples_per_block; /* Maximum number of samples per block */
    uint8_t preset; /* Encoding parameter presets */
};

/* Encoder config */
struct SRLAEncoderConfig {
    uint32_t max_num_channels; /* Maximum number of channels */
    uint32_t min_num_samples_per_block; /* Lower limit of number of samples per block */
    uint32_t max_num_samples_per_block; /* Maximum number of samples per block */
    uint32_t max_num_parameters; /* Maximum number of parameters */
};

/* Encoder handle */
struct SRLAEncoder;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Header encoding */
SRLAApiResult SRLAEncoder_EncodeHeader(
    const struct SRLAHeader *header, uint8_t *data, uint32_t data_size);

/* Calculate the work size required to create the encoder handle */
int32_t SRLAEncoder_CalculateWorkSize(const struct SRLAEncoderConfig *config);

/* Create encoder handle */
struct SRLAEncoder *SRLAEncoder_Create(const struct SRLAEncoderConfig *config, void *work, int32_t work_size);

/* Destroy the encoder handle */
void SRLAEncoder_Destroy(struct SRLAEncoder *encoder);

/* Set encoding parameters */
SRLAApiResult SRLAEncoder_SetEncodeParameter(
    struct SRLAEncoder *encoder, const struct SRLAEncodeParameter *parameter);

/* Calculate the size of a single data block */
SRLAApiResult SRLAEncoder_ComputeBlockSize(
    struct SRLAEncoder *encoder, const int32_t *const *input, uint32_t num_samples,
    uint32_t *output_size);

/* Single data block encoding */
SRLAApiResult SRLAEncoder_EncodeBlock(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size);

/* Encoding including optimal block division search */
SRLAApiResult SRLAEncoder_EncodeOptimalPartitionedBlock(
    struct SRLAEncoder *encoder,
    const int32_t *const *input, uint32_t num_samples,
    uint8_t *data, uint32_t data_size, uint32_t *output_size);

/* Encode the entire file, including the header */
SRLAApiResult SRLAEncoder_EncodeWhole(
    struct SRLAEncoder *encoder,
    const int32_t *const *input, uint32_t num_samples,
    uint8_t *data, uint32_t data_size, uint32_t *output_size, uint8_t variable_block);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* SRLA_ENCODER_H_INCLUDED */
