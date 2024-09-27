#include "srla_decoder.h"

#include <stdlib.h>
#include <string.h>
#include "srla_lpc_synthesize.h"
#include "srla_internal.h"
#include "srla_utility.h"
#include "srla_coder.h"
#include "byte_array.h"
#include "bit_stream.h"
#include "static_huffman.h"
#include "lpc.h"

/* Internal state flag */
#define SRLADECODER_STATUS_FLAG_ALLOCED_BY_OWN  (1 << 0)  /* self-allocated space */
#define SRLADECODER_STATUS_FLAG_SET_HEADER      (1 << 1)  /* Header set */
#define SRLADECODER_STATUS_FLAG_CHECKSUM_CHECK  (1 << 2)  /* Check the checksum */

/* Internal state flag manipulation macro */
#define SRLADECODER_SET_STATUS_FLAG(decoder, flag)    ((decoder->status_flags) |= (flag))
#define SRLADECODER_CLEAR_STATUS_FLAG(decoder, flag)  ((decoder->status_flags) &= ~(flag))
#define SRLADECODER_GET_STATUS_FLAG(decoder, flag)    ((decoder->status_flags) & (flag))

/* Decoder handle */
struct SRLADecoder {
    struct SRLAHeader header; /* Header */
    uint32_t max_num_channels; /* Maximum number of channels that can be decoded */
    uint32_t max_num_parameters; /* Maximum number of parameters */
    struct SRLAPreemphasisFilter **de_emphasis; /* De-emphasis filter */
    int32_t **params_int; /* LPC coefficients for each channel (int) */
    uint32_t *rshifts; /* Right shift amount of LPC coefficients for each channel */
    uint32_t *coef_order; /* LPC coefficient order for each channel */
    const struct StaticHuffmanTree* param_tree; /* Huffman tree of coefficients */
    const struct StaticHuffmanTree* sum_param_tree; /* Huffman tree of summed coefficients */
    const struct SRLAParameterPreset *parameter_preset; /* Parameter presets */
    uint8_t status_flags; /* Internal state flag */
    void *work; /* Work area start pointer */
};

/* Raw data block decode */
static SRLAApiResult SRLADecoder_DecodeRawData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size);
/* Compressed data block decode */
static SRLAApiResult SRLADecoder_DecodeCompressData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size);
/* Silence data block decode */
static SRLAApiResult SRLADecoder_DecodeSilentData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size);

/* Header decode */
SRLAApiResult SRLADecoder_DecodeHeader(
        const uint8_t *data, uint32_t data_size, struct SRLAHeader *header)
{
    const uint8_t *data_pos;
    uint32_t u32buf;
    uint16_t u16buf;
    uint8_t  u8buf;
    struct SRLAHeader tmp_header;

    /* Argument check */
    if ((data == NULL) || (header == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Data size is insufficient */
    if (data_size < SRLA_HEADER_SIZE) {
        return SRLA_APIRESULT_INSUFFICIENT_DATA;
    }

    /* Set the read pointer */
    data_pos = data;

    /* signature */
    {
        uint8_t buf[4];
        ByteArray_GetUint8(data_pos, &buf[0]);
        ByteArray_GetUint8(data_pos, &buf[1]);
        ByteArray_GetUint8(data_pos, &buf[2]);
        ByteArray_GetUint8(data_pos, &buf[3]);
        if ((buf[0] != '1') || (buf[1] != '2')
                || (buf[2] != '4') || (buf[3] != '9')) {
            return SRLA_APIRESULT_INVALID_FORMAT;
        }
    }

    /* If the signature check passes, read to completion without error */

    /* Format version */
    ByteArray_GetUint32BE(data_pos, &u32buf);
    tmp_header.format_version = u32buf;
    /* Encoder version */
    ByteArray_GetUint32BE(data_pos, &u32buf);
    tmp_header.codec_version = u32buf;
    /* Number of channels */
    ByteArray_GetUint16BE(data_pos, &u16buf);
    tmp_header.num_channels = u16buf;
    /* Number of samples */
    ByteArray_GetUint32BE(data_pos, &u32buf);
    tmp_header.num_samples = u32buf;
    /* Sampling rate */
    ByteArray_GetUint32BE(data_pos, &u32buf);
    tmp_header.sampling_rate = u32buf;
    /* Number of bits per sample */
    ByteArray_GetUint16BE(data_pos, &u16buf);
    tmp_header.bits_per_sample = u16buf;
    /* Maximum number of samples per block */
    ByteArray_GetUint32BE(data_pos, &u32buf);
    tmp_header.max_num_samples_per_block = u32buf;
    /* Parameter presets */
    ByteArray_GetUint8(data_pos, &u8buf);
    tmp_header.preset = u8buf;

    /* Check header size */
    SRLA_ASSERT((data_pos - data) == SRLA_HEADER_SIZE);

    /* Successful completion */
    (*header) = tmp_header;
    return SRLA_APIRESULT_OK;
}

/* Check header format */
static SRLAError SRLADecoder_CheckHeaderFormat(const struct SRLAHeader *header)
{
    /* This is an internal module, so if NULL is passed, it will be dropped */
    SRLA_ASSERT(header != NULL);

    /* Format version */
    /* Supplementary Note: For now, if there is a mismatch, an error occurs unconditionally. */
    if (header->format_version != SRLA_FORMAT_VERSION) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Codec version */
    /* Supplementary Note: For now, if there is a mismatch, an error occurs unconditionally. */
    if (header->codec_version != SRLA_CODEC_VERSION) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Number of channels */
    if (header->num_channels == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Number of samples */
    if (header->num_samples == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Sampling rate */
    if (header->sampling_rate == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* bit depth */
    if (header->bits_per_sample == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Maximum number of samples per block */
    if (header->max_num_samples_per_block == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    /* Parameter presets */
    if (header->preset >= SRLA_NUM_PARAMETER_PRESETS) {
        return SRLA_ERROR_INVALID_FORMAT;
    }

    return SRLA_ERROR_OK;
}

/* Calculate the work size required to create a decoder handle */
int32_t SRLADecoder_CalculateWorkSize(const struct SRLADecoderConfig *config)
{
    int32_t work_size;

    /* Argument check */
    if (config == NULL) {
        return -1;
    }

    /* Config check */
    if ((config->max_num_channels == 0)
            || (config->max_num_parameters == 0)) {
        return -1;
    }

    /* Structure size (+ memory alignment) */
    work_size = sizeof(struct SRLADecoder) + SRLA_MEMORY_ALIGNMENT;
    /* De-emphasis filter size */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(struct SRLAPreemphasisFilter, config->max_num_channels, SRLA_NUM_PREEMPHASIS_FILTERS);
    /* Parameter area */
    /* LPC coefficients (int) */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, config->max_num_channels, config->max_num_parameters);
    /* Right shift amount of LPC coefficients for each channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * config->max_num_channels);
    /* LPC coefficient order for each channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * config->max_num_channels);

    return work_size;
}

/* Create a decoder handle */
struct SRLADecoder *SRLADecoder_Create(const struct SRLADecoderConfig *config, void *work, int32_t work_size)
{
    uint32_t ch, l;
    struct SRLADecoder *decoder;
    uint8_t *work_ptr;
    uint8_t tmp_alloc_by_own = 0;

    /* In case of self-allocation of area */
    if ((work == NULL) && (work_size == 0)) {
        if ((work_size = SRLADecoder_CalculateWorkSize(config)) < 0) {
            return NULL;
        }
        work = malloc((uint32_t)work_size);
        tmp_alloc_by_own = 1;
    }

    /* Argument check */
    if ((config == NULL) || (work == NULL)
            || (work_size < SRLADecoder_CalculateWorkSize(config))) {
        return NULL;
    }

    /* Config check */
    if ((config->max_num_channels == 0)
            || (config->max_num_parameters == 0)) {
        return NULL;
    }

    /* Get the work area start pointer */
    work_ptr = (uint8_t *)work;

    /* Allocate structure area */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    decoder = (struct SRLADecoder *)work_ptr;
    work_ptr += sizeof(struct SRLADecoder);

    /* Structure member set */
    decoder->work = work;
    decoder->max_num_channels = config->max_num_channels;
    decoder->max_num_parameters = config->max_num_parameters;
    decoder->status_flags = 0;  /* Clear state */
    if (tmp_alloc_by_own == 1) {
        SRLADECODER_SET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_ALLOCED_BY_OWN);
    }
    if (config->check_checksum == 1) {
        SRLADECODER_SET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_CHECKSUM_CHECK);
    }

    /* Create a de-emphasis filter */
    SRLA_ALLOCATE_2DIMARRAY(decoder->de_emphasis,
            work_ptr, struct SRLAPreemphasisFilter, config->max_num_channels, SRLA_NUM_PREEMPHASIS_FILTERS);

    /* Allocate buffer space and align all pointers */
    /* LPC coefficients (int) */
    SRLA_ALLOCATE_2DIMARRAY(decoder->params_int,
            work_ptr, int32_t, config->max_num_channels, config->max_num_parameters);
    /* Right shift amount of LPC coefficients for each layer */
    work_ptr = (uint8_t*)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    decoder->rshifts = (uint32_t*)work_ptr;
    work_ptr += config->max_num_channels * sizeof(uint32_t);
    /* LPC coefficient order for each layer */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    decoder->coef_order = (uint32_t *)work_ptr;
    work_ptr += config->max_num_channels * sizeof(uint32_t);

    /* Buffer overrun check */
    /* Supplementary Note: There is a possibility that the memory has already been corrupted, so if the check fails, the program will be dropped. */
    SRLA_ASSERT((work_ptr - (uint8_t *)work) <= work_size);

    /* Pre-emphasis filter initialization */
    for (ch = 0; ch < config->max_num_channels; ch++) {
        for (l = 0; l < SRLA_NUM_PREEMPHASIS_FILTERS; l++) {
            SRLAPreemphasisFilter_Initialize(&decoder->de_emphasis[ch][l]);
        }
    }

    /* Create Huffman tree */
    decoder->param_tree = SRLA_GetParameterHuffmanTree();
    decoder->sum_param_tree = SRLA_GetSumParameterHuffmanTree();

    return decoder;
}

/* Destroy the decoder handle */
void SRLADecoder_Destroy(struct SRLADecoder *decoder)
{
    if (decoder != NULL) {
        if (SRLADECODER_GET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_ALLOCED_BY_OWN)) {
            free(decoder->work);
        }
    }
}

/* Set the header in the decoder */
SRLAApiResult SRLADecoder_SetHeader(
        struct SRLADecoder *decoder, const struct SRLAHeader *header)
{
    /* Argument check */
    if ((decoder == NULL) || (header == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check the validity of the header */
    if (SRLADecoder_CheckHeaderFormat(header) != SRLA_ERROR_OK) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }

    /* Check if decoder capacity is exceeded */
    if (decoder->max_num_channels < header->num_channels) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Check maximum number of layers/parameters */
    {
        const struct SRLAParameterPreset* preset = &g_srla_parameter_preset[header->preset];
        if (decoder->max_num_parameters < preset->max_num_parameters) {
            return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
        }
    }

    /* Get encoding preset */
    SRLA_ASSERT(header->preset < SRLA_NUM_PARAMETER_PRESETS);
    decoder->parameter_preset = &g_srla_parameter_preset[header->preset];

    /* Header set */
    decoder->header = (*header);
    SRLADECODER_SET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_SET_HEADER);

    return SRLA_APIRESULT_OK;
}

/* Raw data block decode */
static SRLAApiResult SRLADecoder_DecodeRawData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size)
{
    uint32_t ch, smpl;
    const struct SRLAHeader *header;
    const uint8_t *read_ptr;

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(decoder != NULL);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(data_size > 0);
    SRLA_ASSERT(buffer != NULL);
    SRLA_ASSERT(buffer[0] != NULL);
    SRLA_ASSERT(num_decode_samples > 0);
    SRLA_ASSERT(decode_size != NULL);

    /* Get header */
    header = &(decoder->header);

    /* Assert when there are not enough channels */
    SRLA_ASSERT(num_channels >= header->num_channels);

    /* Check data size */
    if (data_size < (header->bits_per_sample * num_decode_samples * header->num_channels) / 8) {
        return SRLA_APIRESULT_INSUFFICIENT_DATA;
    }

    /* Get raw data with channel interleaving */
    read_ptr = data;
    switch (header->bits_per_sample) {
    case 8:
        for (smpl = 0; smpl < num_decode_samples; smpl++) {
            for (ch = 0; ch < header->num_channels; ch++) {
                uint8_t buf;
                ByteArray_GetUint8(read_ptr, &buf);
                buffer[ch][smpl] = SRLAUTILITY_UINT32_TO_SINT32(buf);
                SRLA_ASSERT((uint32_t)(read_ptr - data) <= data_size);
            }
        }
        break;
    case 16:
        for (smpl = 0; smpl < num_decode_samples; smpl++) {
            for (ch = 0; ch < header->num_channels; ch++) {
                uint16_t buf;
                ByteArray_GetUint16BE(read_ptr, &buf);
                buffer[ch][smpl] = SRLAUTILITY_UINT32_TO_SINT32(buf);
                SRLA_ASSERT((uint32_t)(read_ptr - data) <= data_size);
            }
        }
        break;
    case 24:
        for (smpl = 0; smpl < num_decode_samples; smpl++) {
            for (ch = 0; ch < header->num_channels; ch++) {
                uint32_t buf;
                ByteArray_GetUint24BE(read_ptr, &buf);
                buffer[ch][smpl] = SRLAUTILITY_UINT32_TO_SINT32(buf);
                SRLA_ASSERT((uint32_t)(read_ptr - data) <= data_size);
            }
        }
        break;
    default: SRLA_ASSERT(0);
    }

    /* Get read size */
    (*decode_size) = (uint32_t)(read_ptr - data);

    return SRLA_APIRESULT_OK;
}

/* Compressed data block decode */
static SRLAApiResult SRLADecoder_DecodeCompressData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size)
{
    uint32_t ch;
    int32_t l;
    struct BitStream reader;
    const struct SRLAHeader *header;
    SRLAChannelProcessMethod ch_process_method;

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(decoder != NULL);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(data_size > 0);
    SRLA_ASSERT(buffer != NULL);
    SRLA_ASSERT(buffer[0] != NULL);
    SRLA_ASSERT(num_decode_samples > 0);
    SRLA_ASSERT(decode_size != NULL);

    /* Get header */
    header = &(decoder->header);

    /* Assert when there are not enough channels */
    SRLA_ASSERT(num_channels >= header->num_channels);

    /* Create a bit reader */
    BitReader_Open(&reader, (uint8_t *)data, data_size);

    /* Get multi-channel processing method */
    BitReader_GetBits(&reader, (uint32_t *)&ch_process_method, 2);

    /* Decode parameters */
    /* Pre-emphasis */
    for (ch = 0; ch < num_channels; ch++) {
        uint32_t uval;
        int32_t head;
        /* Pre-emphasis initial value (common to all) */
        BitReader_GetBits(&reader, &uval, header->bits_per_sample + 1U);
        head = SRLAUTILITY_UINT32_TO_SINT32(uval);
        for (l = 0; l < SRLA_NUM_PREEMPHASIS_FILTERS; l++) {
            decoder->de_emphasis[ch][l].prev = head;
        }
        /* Pre-emphasis coefficient */
        for (l = 0; l < SRLA_NUM_PREEMPHASIS_FILTERS; l++) {
            BitReader_GetBits(&reader, &uval, SRLA_PREEMPHASIS_COEF_SHIFT + 1);
            decoder->de_emphasis[ch][l].coef = SRLAUTILITY_UINT32_TO_SINT32(uval);
        }
    }
    /* LPC coefficient order/LPC coefficient right shift amount/LPC coefficient */
    for (ch = 0; ch < num_channels; ch++) {
        uint32_t i, uval, use_sum_coef;
        /* LPC coefficient order */
        BitReader_GetBits(&reader, &decoder->coef_order[ch], SRLA_LPC_COEFFICIENT_ORDER_BITWIDTH);
        decoder->coef_order[ch] += 1; /* Since it has been encoded with -1, return it */
        /* Right shift amount of LPC coefficients in each layer */
        BitReader_GetBits(&reader, &decoder->rshifts[ch], SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH);
        /* LPC coefficients */
        BitReader_GetBits(&reader, &use_sum_coef, 1);
        /* Check whether sum is taken and encoded */
        if (!use_sum_coef) {
            for (i = 0; i < decoder->coef_order[ch]; i++) {
                uval = StaticHuffman_GetCode(decoder->param_tree, &reader);
                decoder->params_int[ch][i] = SRLAUTILITY_UINT32_TO_SINT32(uval);
            }
        } else {
            uval = StaticHuffman_GetCode(decoder->param_tree, &reader);
            decoder->params_int[ch][0] = SRLAUTILITY_UINT32_TO_SINT32(uval);
            for (i = 1; i < decoder->coef_order[ch]; i++) {
                uval = StaticHuffman_GetCode(decoder->sum_param_tree, &reader);
                decoder->params_int[ch][i] = SRLAUTILITY_UINT32_TO_SINT32(uval);
                /* Take the difference and put it back */
                decoder->params_int[ch][i] -= decoder->params_int[ch][i - 1];
            }
        }
    }

    /* Residual decoding */
    for (ch = 0; ch < header->num_channels; ch++) {
        SRLACoder_Decode(&reader, buffer[ch], num_decode_samples);
    }

    /* Align to byte boundary */
    BitStream_Flush(&reader);

    /* Get read size */
    BitStream_Tell(&reader, (int32_t *)decode_size);

    /* Discard bit writer */
    BitStream_Close(&reader);

    /* Composite processing for each channel */
    for (ch = 0; ch < header->num_channels; ch++) {
        /* LPC synthesis */
        SRLALPC_Synthesize(buffer[ch],
            num_decode_samples, decoder->params_int[ch], decoder->coef_order[ch], decoder->rshifts[ch]);
        /* De-emphasis */
        SRLAPreemphasisFilter_MultiStageDeemphasis(
            decoder->de_emphasis[ch], SRLA_NUM_PREEMPHASIS_FILTERS, buffer[ch], num_decode_samples);
    }

    /* Multi-channel processing */
    switch (ch_process_method) {
    case SRLA_CH_PROCESS_METHOD_NONE:
        break;
    case SRLA_CH_PROCESS_METHOD_MS:
        SRLA_ASSERT(header->num_channels >= 2);
        SRLAUtility_MStoLRConversion(buffer, num_decode_samples);
        break;
    case SRLA_CH_PROCESS_METHOD_LS:
        SRLA_ASSERT(header->num_channels >= 2);
        SRLAUtility_LStoLRConversion(buffer, num_decode_samples);
        break;
    case SRLA_CH_PROCESS_METHOD_SR:
        SRLA_ASSERT(header->num_channels >= 2);
        SRLAUtility_SRtoLRConversion(buffer, num_decode_samples);
        break;
    default:
        SRLA_ASSERT(0);
    }

    /* Successful completion */
    return SRLA_APIRESULT_OK;
}

/* Silence data block decode */
static SRLAApiResult SRLADecoder_DecodeSilentData(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t num_channels, uint32_t num_decode_samples,
        uint32_t *decode_size)
{
    uint32_t ch;
    const struct SRLAHeader *header;

    SRLAUTILITY_UNUSED_ARGUMENT(data_size);

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(decoder != NULL);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(buffer != NULL);
    SRLA_ASSERT(buffer[0] != NULL);
    SRLA_ASSERT(num_decode_samples > 0);
    SRLA_ASSERT(decode_size != NULL);

    /* Get header */
    header = &(decoder->header);

    /* Assert when there are not enough channels */
    SRLA_ASSERT(num_channels >= header->num_channels);

    /* Fill everything with silence */
    for (ch = 0; ch < header->num_channels; ch++) {
        memset(buffer[ch], 0, sizeof(int32_t) * num_decode_samples);
    }

    (*decode_size) = 0;
    return SRLA_APIRESULT_OK;
}

/* Single data block decode */
SRLAApiResult SRLADecoder_DecodeBlock(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t buffer_num_channels, uint32_t buffer_num_samples,
        uint32_t *decode_size, uint32_t *num_decode_samples)
{
    uint8_t buf8;
    uint16_t buf16;
    uint32_t buf32;
    uint16_t num_block_samples;
    uint32_t block_header_size, block_data_size;
    SRLAApiResult ret;
    SRLABlockDataType block_type;
    const struct SRLAHeader *header;
    const uint8_t *read_ptr;

    /* Argument check */
    if ((decoder == NULL) || (data == NULL)
            || (buffer == NULL) || (decode_size == NULL)
            || (num_decode_samples == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Header not set yet */
    if (!SRLADECODER_GET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_SET_HEADER)) {
        return SRLA_APIRESULT_PARAMETER_NOT_SET;
    }

    /* Get header */
    header = &(decoder->header);

    /* Check number of buffer channels */
    if (buffer_num_channels < header->num_channels) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Block header decode */
    read_ptr = data;

    /* Synchronous code */
    ByteArray_GetUint16BE(read_ptr, &buf16);
    /* Synchronization code mismatch */
    if (buf16 != SRLA_BLOCK_SYNC_CODE) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* Block size */
    ByteArray_GetUint32BE(read_ptr, &buf32);
    SRLA_ASSERT(buf32 > 0);
    /* Data size is insufficient */
    if ((buf32 + 6) > data_size) {
        return SRLA_APIRESULT_INSUFFICIENT_DATA;
    }
    /* block checksum */
    ByteArray_GetUint16BE(read_ptr, &buf16);
    /* If checking, calculate the checksum and confirm that it matches the obtained value */
    if (SRLADECODER_GET_STATUS_FLAG(decoder, SRLADECODER_STATUS_FLAG_CHECKSUM_CHECK)) {
        /* -2 to exclude the checksum area */
        uint16_t checksum = SRLAUtility_CalculateFletcher16CheckSum(read_ptr, buf32 - 2);
        if (checksum != buf16) {
            return SRLA_APIRESULT_DETECT_DATA_CORRUPTION;
        }
    }
    /* Block data type */
    ByteArray_GetUint8(read_ptr, &buf8);
    block_type = (SRLABlockDataType)buf8;
    /* Number of samples per block channel */
    ByteArray_GetUint16BE(read_ptr, &num_block_samples);
    if (num_block_samples > buffer_num_samples) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }
    /* Block header size */
    block_header_size = (uint32_t)(read_ptr - data);

    /* Decode the data part */
    switch (block_type) {
    case SRLA_BLOCK_DATA_TYPE_RAWDATA:
        ret = SRLADecoder_DecodeRawData(decoder,
                read_ptr, data_size - block_header_size, buffer, header->num_channels, num_block_samples, &block_data_size);
        break;
    case SRLA_BLOCK_DATA_TYPE_COMPRESSDATA:
        ret = SRLADecoder_DecodeCompressData(decoder,
                read_ptr, data_size - block_header_size, buffer, header->num_channels, num_block_samples, &block_data_size);
        break;
    case SRLA_BLOCK_DATA_TYPE_SILENT:
        ret = SRLADecoder_DecodeSilentData(decoder,
                read_ptr, data_size - block_header_size, buffer, header->num_channels, num_block_samples, &block_data_size);
        break;
    default:
        return SRLA_APIRESULT_INVALID_FORMAT;
    }

    /* Data decoding failed */
    if (ret != SRLA_APIRESULT_OK) {
        return ret;
    }

    /* Decode size */
    (*decode_size) = block_header_size + block_data_size;

    /* Number of decoded samples */
    (*num_decode_samples) = num_block_samples;

    /* Decode successful */
    return SRLA_APIRESULT_OK;
}

/* Decode all blocks including header */
SRLAApiResult SRLADecoder_DecodeWhole(
        struct SRLADecoder *decoder,
        const uint8_t *data, uint32_t data_size,
        int32_t **buffer, uint32_t buffer_num_channels, uint32_t buffer_num_samples)
{
    SRLAApiResult ret;
    uint32_t progress, ch, read_offset, read_block_size, num_decode_samples;
    const uint8_t *read_pos;
    int32_t *buffer_ptr[SRLA_MAX_NUM_CHANNELS];
    struct SRLAHeader tmp_header;
    const struct SRLAHeader *header;

    /* Argument check */
    if ((decoder == NULL) || (data == NULL) || (buffer == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Decode the header and set it to the decoder */
    if ((ret = SRLADecoder_DecodeHeader(data, data_size, &tmp_header))
            != SRLA_APIRESULT_OK) {
        return ret;
    }
    if ((ret = SRLADecoder_SetHeader(decoder, &tmp_header))
            != SRLA_APIRESULT_OK) {
        return ret;
    }
    header = &(decoder->header);

    /* Check buffer size */
    if ((buffer_num_channels < header->num_channels)
            || (buffer_num_samples < header->num_samples)) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    progress = 0;
    read_offset = SRLA_HEADER_SIZE;
    read_pos = data + SRLA_HEADER_SIZE;
    while ((progress < header->num_samples) && (read_offset < data_size)) {
        /* Set sample writing position */
        for (ch = 0; ch < header->num_channels; ch++) {
            buffer_ptr[ch] = &buffer[ch][progress];
        }
        /* Block decoding */
        if ((ret = SRLADecoder_DecodeBlock(decoder,
                        read_pos, data_size - read_offset,
                        buffer_ptr, buffer_num_channels, buffer_num_samples - progress,
                        &read_block_size, &num_decode_samples)) != SRLA_APIRESULT_OK) {
            return ret;
        }
        /* Progress update */
        read_pos    += read_block_size;
        read_offset += read_block_size;
        progress    += num_decode_samples;
        SRLA_ASSERT(progress <= buffer_num_samples);
        SRLA_ASSERT(read_offset <= data_size);
    }

    /* Successful completion */
    return SRLA_APIRESULT_OK;
}
