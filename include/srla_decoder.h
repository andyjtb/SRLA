#ifndef SRLA_DECODER_H_INCLUDED
#define SRLA_DECODER_H_INCLUDED

#include "srla.h"
#include "srla_stdint.h"

/* Decoder config */
struct SRLADecoderConfig {
    uint32_t max_num_channels; /* Maximum number of channels */
    uint32_t max_num_parameters; /* Maximum number of parameters */
    uint8_t check_checksum; /* Do you want to check for data corruption using checksum? 1: ON Other: OFF */
};

/* Decoder handle */
struct SRLADecoder;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Header decode */
SRLAApiResult SRLADecoder_DecodeHeader(
        const uint8_t *data, uint32_t data_size, struct SRLAHeader *header);

/* Calculate the work size required to create a decoder handle */
int32_t SRLADecoder_CalculateWorkSize(const struct SRLADecoderConfig *condig);

/* Create a decoder handle */
struct SRLADecoder* SRLADecoder_Create(const struct SRLADecoderConfig *condig, void *work, int32_t work_size);

/* Destroy the decoder handle */
void SRLADecoder_Destroy(struct SRLADecoder *decoder);

/* Set the header in the decoder */
SRLAApiResult SRLADecoder_SetHeader(
        struct SRLADecoder *decoder, const struct SRLAHeader *header);

/* Single data block decode */
SRLAApiResult SRLADecoder_DecodeBlock(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t buffer_num_channels, uint32_t buffer_num_samples,
        uint32_t *decode_size, uint32_t *num_decode_samples);

/* Decode all blocks including header */
SRLAApiResult SRLADecoder_DecodeWhole(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t buffer_num_channels, uint32_t buffer_num_samples);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* SRLA_DECODER_H_INCLUDED */
