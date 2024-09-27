#include "srla_encoder.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "srla_lpc_predict.h"
#include "srla_internal.h"
#include "srla_utility.h"
#include "byte_array.h"
#include "bit_stream.h"
#include "lpc.h"
#include "static_huffman.h"
#include "srla_coder.h"

/* Huge weights when using Dijkstra algorithm */
#define SRLAENCODER_DIJKSTRA_BIGWEIGHT (double)(1UL << 24)

/* Calculate the number of nodes required for block search */
#define SRLAENCODER_CALCULATE_NUM_NODES(num_samples, delta_num_samples) ((SRLAUTILITY_ROUNDUP(num_samples, delta_num_samples) / (delta_num_samples)) + 1)

/* Encoder handle */
struct SRLAEncoder {
    struct SRLAHeader header; /* Header */
    struct SRLACoder *coder; /* encoding handle */
    uint32_t max_num_channels; /* Number of buffer channels */
    uint32_t max_num_samples_per_block; /* Number of buffer samples */
    uint32_t min_num_samples_per_block; /* Minimum number of block samples */
    uint32_t lb_num_samples_per_block; /* Lower limit of number of block samples */
    uint32_t max_num_parameters; /* Maximum number of parameters */
    uint8_t set_parameter; /* Are the parameters set? */
    struct LPCCalculator *lpcc; /* LPC calculation handle */
    struct SRLAPreemphasisFilter **pre_emphasis; /* Pre-emphasis filter */
    struct SRLAOptimalBlockPartitionCalculator *obpc; /* Optimal block division calculation handle */
    double **params_double; /* LPC coefficients for each channel (double) */
    int32_t **params_int; /* LPC coefficients for each channel (int) */
    uint32_t *coef_rshifts; /* Right shift amount of LPC coefficients for each channel */
    uint32_t *coef_order; /* LPC coefficient order for each channel */
    uint32_t *use_sum_coef; /* Is the LPC coefficient summed and encoded for each channel? */
    int32_t **buffer_int; /* signal buffer(int) */
    int32_t **residual; /* Residual signal */
    int32_t **ms_buffer_int; /* MS signal buffer */
    int32_t **ms_residual; /* MS residual signal */
    struct SRLAPreemphasisFilter **ms_pre_emphasis; /* Pre-emphasis filter (MS) */
    uint32_t *ms_coef_rshifts; /* LPC coefficient right shift amount (MS) */
    uint32_t *ms_coef_order; /* LPC coefficient order (MS) */
    uint32_t *ms_use_sum_coef; /* Whether LPC coefficients are summed and coded (MS) */
    int32_t **ms_params_int; /* LPC coefficients (MS) */
    double *buffer_double; /* signal buffer(double) */
    double *error_vars; /* Residual variance sequence for each predictor coefficient */
    double **multiple_lpc_coefs; /* Prediction coefficients for each order */
    uint32_t *partitions_buffer; /* Optimal division setting recording area */
    struct StaticHuffmanCodes param_codes; /* Huffman code for parameter encoding */
    struct StaticHuffmanCodes sum_param_codes; /* Huffman code for parameter encoding taken from the sum */
    const struct SRLAParameterPreset *parameter_preset; /* Parameter presets */
    uint8_t alloced_by_own; /* Is the area allocated by itself? */
    void *work; /* Work area start pointer */
};

/* Optimal block division search handle */
struct SRLAOptimalBlockPartitionCalculator {
    uint32_t max_num_nodes; /* Number of nodes */
    double **adjacency_matrix; /* adjacency matrix */
    double *cost; /* Minimum cost */
    uint32_t *path; /* Path route */
    uint8_t *used_flag; /* Usage status flag for each node */
};

/* Convert encoding parameters to headers */
static SRLAError SRLAEncoder_ConvertParameterToHeader(
        const struct SRLAEncodeParameter *parameter, uint32_t num_samples,
        struct SRLAHeader *header);
/* Determine block data type */
static SRLABlockDataType SRLAEncoder_DecideBlockDataType(
        struct SRLAEncoder *encoder, const int32_t *const *input, uint32_t num_samples);

/* Header encoding */
SRLAApiResult SRLAEncoder_EncodeHeader(
        const struct SRLAHeader *header, uint8_t *data, uint32_t data_size)
{
    uint8_t *data_pos;

    /* Argument check */
    if ((header == NULL) || (data == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Output destination buffer size is insufficient */
    if (data_size < SRLA_HEADER_SIZE) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Check for header abnormal values ​​*/
    /* Do as many checks as possible before writing to data (side effects) */
    /* Number of channels */
    if (header->num_channels == 0) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* Number of samples */
    if (header->num_samples == 0) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* Sampling rate */
    if (header->sampling_rate == 0) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* bit depth */
    if (header->bits_per_sample == 0) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* Maximum number of samples per block */
    if (header->max_num_samples_per_block == 0) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }
    /* Parameter presets */
    if (header->preset >= SRLA_NUM_PARAMETER_PRESETS) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }

    /* Write pointer setting */
    data_pos = data;

    /* signature */
    ByteArray_PutUint8(data_pos, '1');
    ByteArray_PutUint8(data_pos, '2');
    ByteArray_PutUint8(data_pos, '4');
    ByteArray_PutUint8(data_pos, '9');
    /*
/* Format version
* Supplement) Ignores the header setting and writes the macro value */
*/
    ByteArray_PutUint32BE(data_pos, SRLA_FORMAT_VERSION);
    /*
/* Codec version
* Supplement) Ignore the header setting and write the macro value */
*/
    ByteArray_PutUint32BE(data_pos, SRLA_CODEC_VERSION);
    /* Number of channels */
    ByteArray_PutUint16BE(data_pos, header->num_channels);
    /* Number of samples */
    ByteArray_PutUint32BE(data_pos, header->num_samples);
    /* Sampling rate */
    ByteArray_PutUint32BE(data_pos, header->sampling_rate);
    /* Number of bits per sample */
    ByteArray_PutUint16BE(data_pos, header->bits_per_sample);
    /* Maximum number of samples per block */
    ByteArray_PutUint32BE(data_pos, header->max_num_samples_per_block);
    /* Parameter presets */
    ByteArray_PutUint8(data_pos, header->preset);

    /* Check header size */
    SRLA_ASSERT((data_pos - data) == SRLA_HEADER_SIZE);

    /* Successful completion */
    return SRLA_APIRESULT_OK;
}

/* Calculate the work size required to create a search handle */
static int32_t SRLAOptimalBlockPartitionCalculator_CalculateWorkSize(
    uint32_t max_num_samples, uint32_t delta_num_samples)
{
    int32_t work_size;
    uint32_t max_num_nodes;

    /* Calculate the maximum number of nodes */
    max_num_nodes = SRLAENCODER_CALCULATE_NUM_NODES(max_num_samples, delta_num_samples);

    /* Structure size */
    work_size = sizeof(struct SRLAOptimalBlockPartitionCalculator) + SRLA_MEMORY_ALIGNMENT;

    /* adjacency matrix */
    work_size += (int32_t)((sizeof(double *) + (sizeof(double) * max_num_nodes) + SRLA_MEMORY_ALIGNMENT) * max_num_nodes);
    /* Cost array */
    work_size += (int32_t)(sizeof(double) * max_num_nodes + SRLA_MEMORY_ALIGNMENT);
    /* Route information */
    work_size += (int32_t)(sizeof(uint32_t) * max_num_nodes + SRLA_MEMORY_ALIGNMENT);
    /* Node used flag */
    work_size += (int32_t)(sizeof(uint8_t) * max_num_nodes + SRLA_MEMORY_ALIGNMENT);

    return work_size;
}

/* Create a search handle */
static struct SRLAOptimalBlockPartitionCalculator *SRLAOptimalBlockPartitionCalculator_Create(
    uint32_t max_num_samples, uint32_t delta_num_samples, void *work, int32_t work_size)
{
    uint32_t tmp_max_num_nodes;
    struct SRLAOptimalBlockPartitionCalculator* obpc;
    uint8_t *work_ptr;

    /* Argument check */
    if ((max_num_samples < delta_num_samples) || (work == NULL)
        || (work_size < SRLAOptimalBlockPartitionCalculator_CalculateWorkSize(max_num_samples, delta_num_samples))) {
        return NULL;
    }

    /* Get the work area start pointer */
    work_ptr = (uint8_t *)work;

    /* Allocate encoder handle area */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    obpc = (struct SRLAOptimalBlockPartitionCalculator *)work_ptr;
    work_ptr += sizeof(struct SRLAOptimalBlockPartitionCalculator);

    /* Calculate the maximum number of nodes */
    tmp_max_num_nodes = SRLAENCODER_CALCULATE_NUM_NODES(max_num_samples, delta_num_samples);
    obpc->max_num_nodes = tmp_max_num_nodes;

    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    obpc->adjacency_matrix = (double **)work_ptr;
    work_ptr += sizeof(struct SRLAOptimalBlockPartitionCalculator);

    /* Allocate area */
    /* adjacency matrix */
    SRLA_ALLOCATE_2DIMARRAY(obpc->adjacency_matrix, work_ptr, double, tmp_max_num_nodes, tmp_max_num_nodes);
    /* Cost array */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    obpc->cost = (double *)work_ptr;
    work_ptr += sizeof(double) * tmp_max_num_nodes;
    /* Route array */
    work_ptr = (uint8_t*)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    obpc->path = (uint32_t*)work_ptr;
    work_ptr += sizeof(uint32_t) * tmp_max_num_nodes;
    /* Node used flag */
    work_ptr = (uint8_t*)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    obpc->used_flag = (uint8_t*)work_ptr;
    work_ptr += sizeof(uint8_t) * tmp_max_num_nodes;

    return obpc;
}

/* Create a search handle */
static void SRLAOptimalBlockPartitionCalculator_Destroy(struct SRLAOptimalBlockPartitionCalculator *obpc)
{
    /* Do nothing in particular */
    SRLAUTILITY_UNUSED_ARGUMENT(obpc);
}

/* Find the shortest path using Dijkstra's algorithm */
static SRLAError SRLAOptimalBlockPartitionCalculator_ApplyDijkstraMethod(
    struct SRLAOptimalBlockPartitionCalculator *obpc,
    uint32_t num_nodes, uint32_t start_node, uint32_t goal_node, double *min_cost)
{
    uint32_t i, target;
    double min;

    /* Argument check */
    if (obpc == NULL || (num_nodes > obpc->max_num_nodes)) {
        return SRLA_ERROR_INVALID_ARGUMENT;
    }

    /* Clear flags and routes, set distance to a huge value */
    for (i = 0; i < obpc->max_num_nodes; i++) {
        obpc->used_flag[i] = 0;
        obpc->path[i] = ~0U;
        obpc->cost[i] = SRLAENCODER_DIJKSTRA_BIGWEIGHT;
    }

    /* Execute Dijkstra's algorithm */
    obpc->cost[start_node] = 0.0;
    target = start_node;
    while (1) {
        /*
/* The node with the smallest distance (weight) from the undetermined node is determined as the minimum distance to that point. */
*/
        min = SRLAENCODER_DIJKSTRA_BIGWEIGHT;
        for (i = 0; i < num_nodes; i++) {
            if ((obpc->used_flag[i] == 0) && (min > obpc->cost[i])) {
                min = obpc->cost[i];
                target = i;
            }
        }

        /* The shortest route is determined */
        if (target == goal_node) {
            break;
        }

        /*
/* Calculate the distance from the currently determined node to the nodes that are directly connected and not yet determined
* via the currently determined node,
* and if it is smaller than the previous distance, correct the distance and route */
*/
        for (i = 0; i < num_nodes; i++) {
            if (obpc->cost[i] > (obpc->adjacency_matrix[target][i] + obpc->cost[target])) {
                obpc->cost[i] = obpc->adjacency_matrix[target][i] + obpc->cost[target];
                obpc->path[i] = target;
            }
        }

        /* Change the currently focused node to confirmed */
        obpc->used_flag[target] = 1;
    }

    /* Set of minimum costs */
    if (min_cost != NULL) {
        (*min_cost) = obpc->cost[goal_node];
    }

    return SRLA_ERROR_OK;
}

/* Finding the best block division */
static SRLAError SRLAEncoder_SearchOptimalBlockPartitions(
    struct SRLAEncoder *encoder,
    const int32_t *const *input, uint32_t num_samples, uint32_t min_num_block_samples,
    uint32_t *optimal_num_partitions, uint32_t *optimal_block_partition)
{
    uint32_t i, j, ch;
    uint32_t num_channels, num_nodes, tmp_optimal_num_partitions, tmp_node;
    struct SRLAOptimalBlockPartitionCalculator *obpc;

    /* Argument check */
    if ((encoder == NULL) || (input == NULL) || (optimal_num_partitions == NULL)
        || (optimal_block_partition == NULL)) {
        return SRLA_ERROR_INVALID_ARGUMENT;
    }

    /* Received by auto variables */
    obpc = encoder->obpc;
    num_channels = encoder->header.num_channels;

    /* Calculate the adjacency matrix dimension (number of nodes) */
    num_nodes = SRLAENCODER_CALCULATE_NUM_NODES(num_samples, min_num_block_samples);

    /* Maximum number of nodes exceeded */
    if (num_nodes > obpc->max_num_nodes) {
        return SRLA_ERROR_INVALID_ARGUMENT;
    }

    /* Fill the adjacency matrix with huge values ​​once */
    for (i = 0; i < num_nodes; i++) {
        for (j = 0; j < num_nodes; j++) {
            obpc->adjacency_matrix[i][j] = SRLAENCODER_DIJKSTRA_BIGWEIGHT;
        }
    }

    /* Set of adjacency matrices */
    /*
/* The (i,j) element contains the cost (code length) when encoding from i * delta_num_samples to j * delta_num_samples. */
*/
    for (i = 0; i < num_nodes; i++) {
        for (j = i + 1; j < num_nodes; j++) {
            double code_length;
            const int32_t *data_ptr[SRLA_MAX_NUM_CHANNELS];
            const uint32_t sample_offset = i * min_num_block_samples;
            uint32_t num_block_samples = (j - i) * min_num_block_samples;

            /* Adjust as it may pop out at the end points */
            num_block_samples = SRLAUTILITY_MIN(num_block_samples, num_samples - sample_offset);

            /* Set data reference position */
            for (ch = 0; ch < num_channels; ch++) {
                data_ptr[ch] = &input[ch][sample_offset];
            }

            {
                /* Encode and measure length */
                uint32_t encode_len;

                if (SRLAEncoder_ComputeBlockSize(encoder,
                    data_ptr, num_block_samples, &encode_len) != SRLA_APIRESULT_OK) {
                    return SRLA_ERROR_NG;
                }

                code_length = encode_len;
                /* TODO: Try with estimated length too */
            }

            /* Set to adjacency matrix */
            obpc->adjacency_matrix[i][j] = code_length;
        }
    }

    /* Run Dijkstra's algorithm */
    if (SRLAOptimalBlockPartitionCalculator_ApplyDijkstraMethod(
            obpc, num_nodes, 0, num_nodes - 1, NULL) != SRLA_ERROR_OK) {
        return SRLA_ERROR_NG;
    }

    /* Interpretation of the results */
    /* Trace backwards from the goal to the starting position and first find the length of the path */
    tmp_optimal_num_partitions = 0;
    tmp_node = num_nodes - 1;
    while (tmp_node != 0) {
        /* The node traversal order should be in ascending order */
        SRLA_ASSERT(tmp_node > obpc->path[tmp_node]);
        tmp_node = obpc->path[tmp_node];
        tmp_optimal_num_partitions++;
    }

    /* Go through again and set the split size information */
    tmp_node = num_nodes - 1;
    for (i = 0; i < tmp_optimal_num_partitions; i++) {
        const uint32_t sample_offset = obpc->path[tmp_node] * min_num_block_samples;
        uint32_t num_block_samples = (tmp_node - obpc->path[tmp_node]) * min_num_block_samples;
        num_block_samples = SRLAUTILITY_MIN(num_block_samples, num_samples - sample_offset);

        /* Set the number of clipped block samples */
        optimal_block_partition[tmp_optimal_num_partitions - i - 1] = num_block_samples;
        tmp_node = obpc->path[tmp_node];
    }

    /* Set the number of divisions */
    (*optimal_num_partitions) = tmp_optimal_num_partitions;

    return SRLA_ERROR_OK;
}

/* Convert encoding parameters to headers */
static SRLAError SRLAEncoder_ConvertParameterToHeader(
        const struct SRLAEncodeParameter *parameter, uint32_t num_samples,
        struct SRLAHeader *header)
{
    struct SRLAHeader tmp_header = { 0, };

    /* Argument check */
    if ((parameter == NULL) || (header == NULL)) {
        return SRLA_ERROR_INVALID_ARGUMENT;
    }

    /* Check parameters */
    if (parameter->num_channels == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    if (parameter->bits_per_sample == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    if (parameter->sampling_rate == 0) {
        return SRLA_ERROR_INVALID_FORMAT;
    }
    if (parameter->preset >= SRLA_NUM_PARAMETER_PRESETS) {
        return SRLA_ERROR_INVALID_FORMAT;
    }

    /* Total number of samples */
    tmp_header.num_samples = num_samples;

    /* Copy corresponding members */
    tmp_header.num_channels = parameter->num_channels;
    tmp_header.sampling_rate = parameter->sampling_rate;
    tmp_header.bits_per_sample = parameter->bits_per_sample;
    tmp_header.preset = parameter->preset;
    tmp_header.max_num_samples_per_block = parameter->max_num_samples_per_block;

    /* Successful completion */
    (*header) = tmp_header;
    return SRLA_ERROR_OK;
}

/* Calculate the work size required to create the encoder handle */
int32_t SRLAEncoder_CalculateWorkSize(const struct SRLAEncoderConfig *config)
{
    int32_t work_size, tmp_work_size;

    /* Argument check */
    if (config == NULL) {
        return -1;
    }

    /* Config check */
    if ((config->max_num_samples_per_block == 0)
            || (config->min_num_samples_per_block == 0)
            || (config->max_num_channels == 0)
            || (config->max_num_parameters == 0)) {
        return -1;
    }

    /* Block size should be greater than the number of parameters */
    if (config->max_num_parameters > config->max_num_samples_per_block) {
        return -1;
    }
    /* Lower block size limit exceeds upper limit */
    if (config->min_num_samples_per_block > config->max_num_samples_per_block) {
        return -1;
    }

    /* Handle body size */
    work_size = sizeof(struct SRLAEncoder) + SRLA_MEMORY_ALIGNMENT;

    /* Size of LPC calculation handle */
    {
        struct LPCCalculatorConfig lpcc_config;
        lpcc_config.max_num_samples = config->max_num_samples_per_block;
        lpcc_config.max_order = config->max_num_parameters;
        if ((tmp_work_size = LPCCalculator_CalculateWorkSize(&lpcc_config)) < 0) {
            return -1;
        }
        work_size += tmp_work_size;
    }

    /* Size of the encoded handle */
    if ((tmp_work_size = SRLACoder_CalculateWorkSize(config->max_num_samples_per_block)) < 0) {
        return -1;
    }
    work_size += tmp_work_size;

    /* Size of optimal split search handle */
    if ((tmp_work_size = SRLAOptimalBlockPartitionCalculator_CalculateWorkSize(
            config->max_num_samples_per_block, config->min_num_samples_per_block)) < 0) {
        return -1;
    }
    work_size += tmp_work_size;

    /* Pre-emphasis filter size */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(struct SRLAPreemphasisFilter, config->max_num_channels, SRLA_NUM_PREEMPHASIS_FILTERS);
    /* Pre-emphasis filter size (MS) */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(struct SRLAPreemphasisFilter, 2, SRLA_NUM_PREEMPHASIS_FILTERS);
    /* Parameter buffer area */
    /* LPC coefficients (int) */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, config->max_num_channels, config->max_num_parameters);
    /* LPC coefficients (MS, int) */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, 2, config->max_num_parameters);
    /* LPC coefficients (double) */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(double, config->max_num_channels, config->max_num_parameters);
    /* Right shift amount of LPC coefficients for each channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * config->max_num_channels);
    /* Right shift amount of LPC coefficients for MS channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * 2);
    /* LPC coefficient order for each channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * config->max_num_channels);
    /* LPC coefficient order for MS channel */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * 2);
    /* Flag indicating whether LPC coefficients for each channel are summed and encoded */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * config->max_num_channels);
    /* Flag indicating whether the LPC coefficients of the MS channel are summed and encoded */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(uint32_t) * 2);
    /* Size of signal processing buffer */
    work_size += (int32_t)(2 * SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, config->max_num_channels, config->max_num_samples_per_block));
    work_size += (int32_t)(config->max_num_samples_per_block * sizeof(double) + SRLA_MEMORY_ALIGNMENT);
    /* Size of the residual signal */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, config->max_num_channels, config->max_num_samples_per_block);
    /* MS signal and residual buffer size */
    work_size += (int32_t)(2 * SRLA_CALCULATE_2DIMARRAY_WORKSIZE(int32_t, 2, config->max_num_samples_per_block));
    /* Size of residual variance region */
    work_size += (int32_t)(SRLA_MEMORY_ALIGNMENT + sizeof(double) * (config->max_num_parameters + 1));
    /* Size of coefficient area */
    work_size += (int32_t)SRLA_CALCULATE_2DIMARRAY_WORKSIZE(double, config->max_num_parameters, config->max_num_parameters);
    /* Size of split recording area */
    work_size += (int32_t)(SRLAENCODER_CALCULATE_NUM_NODES(config->max_num_samples_per_block, config->min_num_samples_per_block) * sizeof(uint32_t) + SRLA_MEMORY_ALIGNMENT);

    return work_size;
}

/* Create encoder handle */
struct SRLAEncoder* SRLAEncoder_Create(const struct SRLAEncoderConfig *config, void *work, int32_t work_size)
{
    uint32_t ch, l;
    struct SRLAEncoder *encoder;
    uint8_t tmp_alloc_by_own = 0;
    uint8_t *work_ptr;

    /* In case of pre-allocation of work area */
    if ((work == NULL) && (work_size == 0)) {
        if ((work_size = SRLAEncoder_CalculateWorkSize(config)) < 0) {
            return NULL;
        }
        work = malloc((uint32_t)work_size);
        tmp_alloc_by_own = 1;
    }

    /* Argument check */
    if ((config == NULL) || (work == NULL)
        || (work_size < SRLAEncoder_CalculateWorkSize(config))) {
        return NULL;
    }

    /* Config check */
    if ((config->max_num_channels == 0)
        || (config->max_num_samples_per_block == 0)
        || (config->max_num_parameters == 0)) {
        return NULL;
    }

    /* Block size should be greater than the number of parameters */
    if (config->max_num_parameters > config->max_num_samples_per_block) {
        return NULL;
    }

    /* Get the work area start pointer */
    work_ptr = (uint8_t *)work;

    /* Allocate encoder handle area */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder = (struct SRLAEncoder *)work_ptr;
    work_ptr += sizeof(struct SRLAEncoder);

    /* Encoder member settings */
    encoder->set_parameter = 0;
    encoder->alloced_by_own = tmp_alloc_by_own;
    encoder->work = work;
    encoder->max_num_channels = config->max_num_channels;
    encoder->max_num_samples_per_block = config->max_num_samples_per_block;
    encoder->lb_num_samples_per_block = config->min_num_samples_per_block;
    encoder->max_num_parameters = config->max_num_parameters;

    /* Create LPC calculation handle */
    {
        int32_t lpcc_size;
        struct LPCCalculatorConfig lpcc_config;
        lpcc_config.max_num_samples = config->max_num_samples_per_block;
        lpcc_config.max_order = config->max_num_parameters;
        lpcc_size = LPCCalculator_CalculateWorkSize(&lpcc_config);
        if ((encoder->lpcc = LPCCalculator_Create(&lpcc_config, work_ptr, lpcc_size)) == NULL) {
            return NULL;
        }
        work_ptr += lpcc_size;
    }

    /* Create an encoding handle */
    {
        const int32_t coder_size = SRLACoder_CalculateWorkSize(config->max_num_samples_per_block);
        if ((encoder->coder = SRLACoder_Create(config->max_num_samples_per_block, work_ptr, coder_size)) == NULL) {
            return NULL;
        }
        work_ptr += coder_size;
    }

    /* Create optimal split search handle */
    {
        const int32_t obpc_size
            = SRLAOptimalBlockPartitionCalculator_CalculateWorkSize(
            config->max_num_samples_per_block, config->min_num_samples_per_block);
        if ((encoder->obpc = SRLAOptimalBlockPartitionCalculator_Create(
                config->max_num_samples_per_block, config->min_num_samples_per_block, work_ptr, obpc_size)) == NULL) {
            return NULL;
        }
        work_ptr += obpc_size;
    }

    /* Create a pre-emphasis filter */
    SRLA_ALLOCATE_2DIMARRAY(encoder->pre_emphasis,
        work_ptr, struct SRLAPreemphasisFilter, config->max_num_channels, SRLA_NUM_PREEMPHASIS_FILTERS);
    /* Create a pre-emphasis filter for MS channels */
    SRLA_ALLOCATE_2DIMARRAY(encoder->ms_pre_emphasis,
        work_ptr, struct SRLAPreemphasisFilter, 2, SRLA_NUM_PREEMPHASIS_FILTERS);

    /* Allocate buffer space and align all pointers */
    /* LPC coefficients (int) */
    SRLA_ALLOCATE_2DIMARRAY(encoder->params_int,
        work_ptr, int32_t, config->max_num_channels, config->max_num_parameters);
    /* LPC coefficients (MS, int) */
    SRLA_ALLOCATE_2DIMARRAY(encoder->ms_params_int,
        work_ptr, int32_t, 2, config->max_num_parameters);
    /* LPC coefficients (double) */
    SRLA_ALLOCATE_2DIMARRAY(encoder->params_double,
        work_ptr, double, config->max_num_channels, config->max_num_parameters);
    /* LPC coefficient right shift amount */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->coef_rshifts = (uint32_t *)work_ptr;
    work_ptr += config->max_num_channels * sizeof(uint32_t);
    /* LPC coefficient right shift amount (MS) */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->ms_coef_rshifts = (uint32_t *)work_ptr;
    work_ptr += 2 * sizeof(uint32_t);
    /* LPC coefficient order */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->coef_order = (uint32_t *)work_ptr;
    work_ptr += config->max_num_channels * sizeof(uint32_t);
    /* LPC coefficient order (MS) */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->ms_coef_order = (uint32_t *)work_ptr;
    work_ptr += 2 * sizeof(uint32_t);
    /* Flag for whether LPC coefficients are summed and coded */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->use_sum_coef = (uint32_t *)work_ptr;
    work_ptr += config->max_num_channels * sizeof(uint32_t);
    /* Flag indicating whether LPC coefficients are summed and coded (MS) */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->ms_use_sum_coef = (uint32_t *)work_ptr;
    work_ptr += 2 * sizeof(uint32_t);

    /* Buffer area for signal processing */
    SRLA_ALLOCATE_2DIMARRAY(encoder->buffer_int,
            work_ptr, int32_t, config->max_num_channels, config->max_num_samples_per_block);
    SRLA_ALLOCATE_2DIMARRAY(encoder->residual,
            work_ptr, int32_t, config->max_num_channels, config->max_num_samples_per_block);

    /* Buffer area for MS signal processing */
    SRLA_ALLOCATE_2DIMARRAY(encoder->ms_buffer_int,
        work_ptr, int32_t, 2, config->max_num_samples_per_block);
    SRLA_ALLOCATE_2DIMARRAY(encoder->ms_residual,
        work_ptr, int32_t, 2, config->max_num_samples_per_block);

    /* Residual variance domain */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->error_vars = (double *)work_ptr;
    work_ptr += (config->max_num_parameters + 1) * sizeof(double);

    /* previous order coefficient */
    SRLA_ALLOCATE_2DIMARRAY(encoder->multiple_lpc_coefs,
        work_ptr, double, config->max_num_parameters, config->max_num_parameters);

    /* double buffer */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->buffer_double = (double *)work_ptr;
    work_ptr += config->max_num_samples_per_block * sizeof(double);

    /* Split setting recording area */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);
    encoder->partitions_buffer = (uint32_t *)work_ptr;
    work_ptr += SRLAENCODER_CALCULATE_NUM_NODES(config->max_num_samples_per_block, config->min_num_samples_per_block) * sizeof(uint32_t);

    /* Buffer overrun check */
    /* Supplementary Note: There is a possibility that the memory has already been corrupted, so if the check fails, the program will be dropped. */
    SRLA_ASSERT((work_ptr - (uint8_t *)work) <= work_size);

    /* Pre-emphasis filter initialization */
    for (ch = 0; ch < config->max_num_channels; ch++) {
        for (l = 0; l < SRLA_NUM_PREEMPHASIS_FILTERS; l++) {
            SRLAPreemphasisFilter_Initialize(&encoder->pre_emphasis[ch][l]);
        }
    }

    /* Create Huffman code */
    StaticHuffman_ConvertTreeToCodes(SRLA_GetParameterHuffmanTree(), &encoder->param_codes);
    StaticHuffman_ConvertTreeToCodes(SRLA_GetSumParameterHuffmanTree(), &encoder->sum_param_codes);

    return encoder;
}

/* Destroy the encoder handle */
void SRLAEncoder_Destroy(struct SRLAEncoder *encoder)
{
    if (encoder != NULL) {
        SRLACoder_Destroy(encoder->coder);
        SRLAOptimalBlockPartitionCalculator_Destroy(encoder->obpc);
        LPCCalculator_Destroy(encoder->lpcc);
        if (encoder->alloced_by_own == 1) {
            free(encoder->work);
        }
    }
}

/* Set encoding parameters */
SRLAApiResult SRLAEncoder_SetEncodeParameter(
        struct SRLAEncoder *encoder, const struct SRLAEncodeParameter *parameter)
{
    struct SRLAHeader tmp_header;

    /* Argument check */
    if ((encoder == NULL) || (parameter == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check whether the parameter settings are correct by converting them to headers */
    /* Enter a dummy value for the total number of samples */
    if (SRLAEncoder_ConvertParameterToHeader(parameter, 0, &tmp_header) != SRLA_ERROR_OK) {
        return SRLA_APIRESULT_INVALID_FORMAT;
    }

    /* Check if the encoder capacity is exceeded */
    if ((encoder->max_num_samples_per_block < parameter->max_num_samples_per_block)
        || (encoder->lb_num_samples_per_block > parameter->min_num_samples_per_block)
        || (encoder->max_num_channels < parameter->num_channels)) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }
    /* Set maximum number of samples per block */
    tmp_header.max_num_samples_per_block = parameter->max_num_samples_per_block;
    /* Record the minimum number of samples per block */
    encoder->min_num_samples_per_block = parameter->min_num_samples_per_block;

    /* Header settings */
    encoder->header = tmp_header;

    /* Get encoding preset */
    SRLA_ASSERT(parameter->preset < SRLA_NUM_PARAMETER_PRESETS);
    encoder->parameter_preset = &g_srla_parameter_preset[parameter->preset];

    /* Set the parameter set flag */
    encoder->set_parameter = 1;

    return SRLA_APIRESULT_OK;
}

/* Determine block data type */
static SRLABlockDataType SRLAEncoder_DecideBlockDataType(
        struct SRLAEncoder *encoder, const int32_t *const *input, uint32_t num_samples)
{
    uint32_t ch, smpl;
    const struct SRLAHeader *header;

    SRLA_ASSERT(encoder != NULL);
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(encoder->set_parameter == 1);

    /* If the order is less than LPC, use raw data */
    if (num_samples <= encoder->parameter_preset->max_num_parameters) {
        return SRLA_BLOCK_DATA_TYPE_RAWDATA;
    }

    header = &encoder->header;

    /* Silence detection */
    for (ch = 0; ch < header->num_channels; ch++) {
        for (smpl = 0; smpl < num_samples; smpl++) {
            if (input[ch][smpl] != 0) {
                goto NOT_SILENCE;
            }
        }
    }
    return SRLA_BLOCK_DATA_TYPE_SILENT;

NOT_SILENCE:
    /* The rest is compressed data */
    return SRLA_BLOCK_DATA_TYPE_COMPRESSDATA;
}

/* Raw data block encoding */
static SRLAApiResult SRLAEncoder_EncodeRawData(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size)
{
    uint32_t ch, smpl;
    const struct SRLAHeader *header;
    uint8_t *data_ptr;

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(encoder != NULL);
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(num_samples > 0);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(data_size > 0);
    SRLA_ASSERT(output_size != NULL);

    header = &(encoder->header);

    /* Check the buffer size of the write destination */
    if (data_size < (header->bits_per_sample * num_samples * header->num_channels) / 8) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Output raw data with channel interleaving */
    data_ptr = data;
    switch (header->bits_per_sample) {
            case 8:
                for (smpl = 0; smpl < num_samples; smpl++) {
                    for (ch = 0; ch < header->num_channels; ch++) {
                        ByteArray_PutUint8(data_ptr, SRLAUTILITY_SINT32_TO_UINT32(input[ch][smpl]));
                        SRLA_ASSERT((uint32_t)(data_ptr - data) < data_size);
                    }
                }
                break;
            case 16:
                for (smpl = 0; smpl < num_samples; smpl++) {
                    for (ch = 0; ch < header->num_channels; ch++) {
                        ByteArray_PutUint16BE(data_ptr, SRLAUTILITY_SINT32_TO_UINT32(input[ch][smpl]));
                        SRLA_ASSERT((uint32_t)(data_ptr - data) < data_size);
                    }
                }
                break;
            case 24:
                for (smpl = 0; smpl < num_samples; smpl++) {
                    for (ch = 0; ch < header->num_channels; ch++) {
                        ByteArray_PutUint24BE(data_ptr, SRLAUTILITY_SINT32_TO_UINT32(input[ch][smpl]));
                        SRLA_ASSERT((uint32_t)(data_ptr - data) < data_size);
                    }
                }
                break;
            default:
                SRLA_ASSERT(0);
    }

    /* Get write size */
    (*output_size) = (uint32_t)(data_ptr - data);

    return SRLA_APIRESULT_OK;
}

/* Average code length of Recursive Golomb-Rice code */
static double SRLAEncoder_CalculateRGRMeanCodeLength(double mean_abs_error, uint32_t bps)
{
    const double intmean = mean_abs_error * (1 << (bps - 1)); /* Average value when integer quantized */
    const double rho = 1.0 / (1.0 + intmean);
    const uint32_t k2 = (uint32_t)SRLAUTILITY_MAX(0, SRLAUtility_Log2(log(0.5127629514) / log(1.0 - rho)));
    const uint32_t k1 = k2 + 1;
    const double k1factor = pow(1.0 - rho, (double)(1 << k1));
    const double k2factor = pow(1.0 - rho, (double)(1 << k2));
    return (1.0 + k1) * (1.0 - k1factor) + (1.0 + k2 + (1.0 / (1.0 - k2factor))) * k1factor;
}

/* Calculate the entropy of the geometric distribution */
static double SRLAEncoder_CalculateGeometricDistributionEntropy(double mean_abs_error, uint32_t bps)
{
    const double intmean = mean_abs_error * (1 << (bps - 1)); /* Average value when integer quantized */
    const double rho = 1.0 / (1.0 + intmean);
    const double invrho = 1.0 - rho;

    if (mean_abs_error < FLT_MIN) {
        return 0.0;
    }

    return -(invrho * SRLAUtility_Log2(invrho) + rho * SRLAUtility_Log2(rho)) / rho;
}

/* Select the optimal LPC order */
static SRLAError SRLAEncoder_SelectBestLPCOrder(
    const struct SRLAHeader *header, SRLAChannelLPCOrderDecisionTactics tactics,
    const double *input, uint32_t num_samples, const double **coefs, const double *error_vars,
    uint32_t max_coef_order, uint32_t *best_coef_order)
{
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(coefs != NULL);
    SRLA_ASSERT(error_vars != NULL);
    SRLA_ASSERT(best_coef_order != NULL);

    switch (tactics) {
    case SRLA_LPC_ORDER_DECISION_TACTICS_MAX_FIXED:
        /* Always choose the highest order */
        (*best_coef_order) = max_coef_order;
        return SRLA_ERROR_OK;
    case SRLA_LPC_ORDER_DECISION_TACTICS_BRUTEFORCE_SEARCH:
        /* exhaustive search */
    {
        double minlen, len, mabse;
        uint32_t i, order, smpl, tmp_best_order = 0;

        minlen = FLT_MAX;
        for (order = 1; order <= max_coef_order; order++) {
            const double *coef = coefs[order - 1];
            mabse = 0.0;
            for (smpl = order; smpl < num_samples; smpl++) {
                double residual = input[smpl];
                for (i = 0; i < order; i++) {
                    residual += coef[i] * input[smpl - i - 1];
                }
                mabse += SRLAUTILITY_ABS(residual);
            }
            /* Size of residual code is doubled to make it a non-negative integer during encoding */
            len = SRLAEncoder_CalculateRGRMeanCodeLength(2.0 * mabse / num_samples, header->bits_per_sample) * num_samples;
            /* Size of coefficient */
            len += SRLA_LPC_COEFFICIENT_BITWIDTH * order;
            if (minlen > len) {
                minlen = len;
                tmp_best_order = order;
            }
        }
        /* Set the result */
        SRLA_ASSERT(tmp_best_order != 0);
        (*best_coef_order) = tmp_best_order;
        return SRLA_ERROR_OK;
    }
    case SRLA_LPC_ORDER_DECISION_TACTICS_BRUTEFORCE_ESTIMATION:
        /* Search for coefficient order that gives minimum estimated code length */
    {
        uint32_t order, tmp_best_order = 0;
        double len, mabse, minlen = FLT_MAX;

        for (order = 1; order <= max_coef_order; order++) {
            /* Estimate the mean absolute value from the residual variance under the assumption of Laplace distribution */
            mabse = 2.0 * sqrt(error_vars[order] / 2.0); /* Doubled to make it a non-negative integer when encoding */
            /* size of residual code */
            len = SRLAEncoder_CalculateGeometricDistributionEntropy(mabse, header->bits_per_sample) * num_samples;
            /* Size of coefficient */
            len += SRLA_LPC_COEFFICIENT_BITWIDTH * order;
            if (minlen > len) {
                minlen = len;
                tmp_best_order = order;
            }
        }

        /* Set the result */
        SRLA_ASSERT(tmp_best_order != 0);
        (*best_coef_order) = tmp_best_order;
        return SRLA_ERROR_OK;
    }
    default:
        SRLA_ASSERT(0);
    }

    return SRLA_ERROR_NG;
}

/* Calculate parameters and code length for one channel */
static SRLAError SRLAEncoder_ComputeCoefficientsPerChannel(
    struct SRLAEncoder *encoder,
    int32_t *buffer_int, double *buffer_double, int32_t *residual_int, uint32_t num_samples,
    struct SRLAPreemphasisFilter *pre_emphasis_filters, uint32_t *coef_order, uint32_t *coef_rshift,
    int32_t *int_coef, uint32_t *use_sum_coef, uint32_t *code_length)
{
    uint32_t smpl, p;
    LPCApiResult ret;
    SRLAError err;
    const struct SRLAHeader *header;
    const struct SRLAParameterPreset *parameter_preset;
    struct SRLAPreemphasisFilter tmp_pre_emphasis_filters[SRLA_NUM_PREEMPHASIS_FILTERS];
    uint32_t tmp_coef_order;
    uint32_t tmp_coef_rshift;
    double *double_coef;
    int32_t tmp_int_coef[SRLA_MAX_COEFFICIENT_ORDER];
    uint32_t tmp_use_sum_coef;
    uint32_t tmp_code_length;

    /* Argument check */
    if ((encoder == NULL) || (buffer_int == NULL) || (buffer_double == NULL) || (residual_int == NULL)
        || (pre_emphasis_filters == NULL) || (coef_order == NULL) || (coef_rshift == NULL) || (int_coef == NULL)
        || (use_sum_coef == NULL) || (code_length == NULL)) {
        return SRLA_ERROR_INVALID_ARGUMENT;
    }

    header = &(encoder->header);
    parameter_preset = encoder->parameter_preset;

    /* Pre-emphasis filters */
    {
        const int32_t head = buffer_int[0];
        struct SRLAPreemphasisFilter filter[SRLA_NUM_PREEMPHASIS_FILTERS] = { 0, };
        SRLAPreemphasisFilter_CalculateMultiStageCoefficients(filter, SRLA_NUM_PREEMPHASIS_FILTERS, buffer_int, num_samples);
        for (p = 0; p < SRLA_NUM_PREEMPHASIS_FILTERS; p++) {
            filter[p].prev = head;
            SRLAPreemphasisFilter_Preemphasis(&filter[p], buffer_int, num_samples);
            tmp_pre_emphasis_filters[p].prev = head;
            tmp_pre_emphasis_filters[p].coef = filter[p].coef;
        }
    }

    /* Convert to a double precision signal (normalized to the range [-1,1]) */
    {
        const double norm_const = pow(2.0, -(int32_t)(header->bits_per_sample - 1));
        for (smpl = 0; smpl < num_samples; smpl++) {
            buffer_double[smpl] = buffer_int[smpl] * norm_const;
        }
    }

    /* Calculate coefficients and error variance up to maximum order */
    if ((ret = LPCCalculator_CalculateMultipleLPCCoefficients(encoder->lpcc,
        buffer_double, num_samples,
        encoder->multiple_lpc_coefs, encoder->error_vars, parameter_preset->max_num_parameters,
        LPC_WINDOWTYPE_WELCH, SRLA_LPC_RIDGE_REGULARIZATION_PARAMETER)) != LPC_APIRESULT_OK) {
        return SRLA_ERROR_NG;
    };

    /* Select degree */
    if ((err = SRLAEncoder_SelectBestLPCOrder(header,
        parameter_preset->lpc_order_tactics,
        buffer_double, num_samples, (const double **)encoder->multiple_lpc_coefs, encoder->error_vars,
        parameter_preset->max_num_parameters, &tmp_coef_order)) != SRLA_ERROR_OK) {
        return err;
    }

    /* Set the best order as a parameter */
    double_coef = encoder->multiple_lpc_coefs[tmp_coef_order - 1];

    /* LPC coefficient calculation using SVR */
    if ((ret = LPCCalculator_CalculateLPCCoefficientsSVR(encoder->lpcc,
        buffer_double, num_samples,
        double_coef, tmp_coef_order, parameter_preset->svr_max_num_iterations,
        LPC_WINDOWTYPE_WELCH, SRLA_LPC_RIDGE_REGULARIZATION_PARAMETER,
        parameter_preset->margin_list, parameter_preset->margin_list_size)) != LPC_APIRESULT_OK) {
        return SRLA_ERROR_NG;
    }

    /* LPC coefficient quantization */
    if ((ret = LPC_QuantizeCoefficients(double_coef, tmp_coef_order,
        SRLA_LPC_COEFFICIENT_BITWIDTH, (1 << SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH),
        tmp_int_coef, &tmp_coef_rshift)) != LPC_APIRESULT_OK) {
        return SRLA_ERROR_NG;
    }

    /* Change the parameter order to increase the index in the convolution operation */
    for (p = 0; p < tmp_coef_order / 2; p++) {
        int32_t tmp = tmp_int_coef[p];
        tmp_int_coef[p] = tmp_int_coef[tmp_coef_order - p - 1];
        tmp_int_coef[tmp_coef_order - p - 1] = tmp;
    }

    /* LPC prediction */
    SRLALPC_Predict(buffer_int,
        num_samples, tmp_int_coef, tmp_coef_order, residual_int, tmp_coef_rshift);

    /* Code length calculation */
    tmp_code_length = 0;

    /* Residual code length */
    tmp_code_length += SRLACoder_ComputeCodeLength(encoder->coder, residual_int, num_samples);

    /* Pre-emphasis filter buffer/coefficients */
    tmp_code_length += header->bits_per_sample + 1;
    for (p = 0; p < SRLA_NUM_PREEMPHASIS_FILTERS; p++) {
        tmp_code_length += SRLA_PREEMPHASIS_COEF_SHIFT + 1;
    }

    /* LPC coefficient order/LPC coefficient right shift amount */
    tmp_code_length += SRLA_LPC_COEFFICIENT_ORDER_BITWIDTH;
    tmp_code_length += SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH;

    /* Flag field indicating whether or not sum was taken */
    tmp_code_length += 1;

    /* LPC coefficients area */
    {
        uint32_t coef_code_length = 0, summed_coef_code_length;

        /* Calculate code length without summing */
        for (p = 0; p < tmp_coef_order; p++) {
            const uint32_t uval = SRLAUTILITY_SINT32_TO_UINT32(tmp_int_coef[p]);
            SRLA_ASSERT(uval < STATICHUFFMAN_MAX_NUM_SYMBOLS);
            coef_code_length += encoder->param_codes.codes[uval].bit_count;
        }

        /* Calculate the code length by taking the sum */
        tmp_use_sum_coef = 1;
        summed_coef_code_length
            = encoder->param_codes.codes[SRLAUTILITY_SINT32_TO_UINT32(tmp_int_coef[0])].bit_count;
        for (p = 1; p < tmp_coef_order; p++) {
            const int32_t summed = tmp_int_coef[p] + tmp_int_coef[p - 1];
            const uint32_t uval = SRLAUTILITY_SINT32_TO_UINT32(summed);
            if (uval >= STATICHUFFMAN_MAX_NUM_SYMBOLS) {
                tmp_use_sum_coef = 0;
                break;
            }
            summed_coef_code_length += encoder->sum_param_codes.codes[uval].bit_count;
            if (summed_coef_code_length >= coef_code_length) {
                tmp_use_sum_coef = 0;
                break;
            }
        }

        tmp_code_length += (tmp_use_sum_coef) ? summed_coef_code_length : coef_code_length;
    }

    /* Output the results */
    memcpy(pre_emphasis_filters, tmp_pre_emphasis_filters,
        sizeof(struct SRLAPreemphasisFilter) * SRLA_NUM_PREEMPHASIS_FILTERS);
    (*coef_order) = tmp_coef_order;
    (*coef_rshift) = tmp_coef_rshift;
    memcpy(int_coef, tmp_int_coef, sizeof(int32_t) * tmp_coef_order);
    (*use_sum_coef) = tmp_use_sum_coef;
    (*code_length) = tmp_code_length;

    return SRLA_ERROR_OK;
}

/* Calculate coefficients for compressed data block */
static SRLAApiResult SRLAEncoder_ComputeCoefficients(
    struct SRLAEncoder *encoder, const int32_t *const *input, uint32_t num_samples,
    SRLAChannelProcessMethod *ch_process_method, uint32_t *output_bits)
{
    uint32_t ch, tmp_output_bits = 0;
    const struct SRLAHeader *header;
    SRLAChannelProcessMethod tmp_ch_process_method = SRLA_CH_PROCESS_METHOD_INVALID;
    uint32_t code_length[SRLA_MAX_NUM_CHANNELS] = { 0, };
    uint32_t ms_code_length[2] = { 0, };

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(encoder != NULL);
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(num_samples > 0);
    SRLA_ASSERT(ch_process_method != NULL);
    SRLA_ASSERT(output_bits != NULL);

    /* Get header */
    header = &(encoder->header);

    /* Copy input to buffer */
    for (ch = 0; ch < header->num_channels; ch++) {
        memcpy(encoder->buffer_int[ch], input[ch], sizeof(int32_t) * num_samples);
        /* If the input is smaller than the buffer size, pad the end with zeros */
        if (num_samples < encoder->max_num_samples_per_block) {
            const uint32_t remain = encoder->max_num_samples_per_block - num_samples;
            memset(&encoder->buffer_int[ch][num_samples], 0, sizeof(int32_t) * remain);
        }
    }

    /* MS signal generation/code length calculation */
    if (header->num_channels >= 2) {
        for (ch = 0; ch < 2; ch++) {
            memcpy(encoder->ms_buffer_int[ch], encoder->buffer_int[ch], sizeof(int32_t) * num_samples);
        }
        SRLAUtility_LRtoMSConversion(encoder->ms_buffer_int, num_samples);
        for (ch = 0; ch < 2; ch++) {
            SRLAError err;
            if ((err = SRLAEncoder_ComputeCoefficientsPerChannel(encoder,
                encoder->ms_buffer_int[ch], encoder->buffer_double, encoder->ms_residual[ch], num_samples,
                encoder->ms_pre_emphasis[ch], &encoder->ms_coef_order[ch], &encoder->ms_coef_rshifts[ch],
                encoder->ms_params_int[ch], &encoder->ms_use_sum_coef[ch], &ms_code_length[ch])) != SRLA_ERROR_OK) {
                return SRLA_APIRESULT_NG;
            }
        }
    }
    /* Calculate parameters and code length for each channel */
    for (ch = 0; ch < header->num_channels; ch++) {
        SRLAError err;
        if ((err = SRLAEncoder_ComputeCoefficientsPerChannel(encoder,
            encoder->buffer_int[ch], encoder->buffer_double, encoder->residual[ch], num_samples,
            encoder->pre_emphasis[ch], &encoder->coef_order[ch], &encoder->coef_rshifts[ch],
            encoder->params_int[ch], &encoder->use_sum_coef[ch], &code_length[ch])) != SRLA_ERROR_OK) {
            return SRLA_APIRESULT_NG;
        }
    }

    /* Select the smallest code length in multi-channel */
    if (header->num_channels == 1) {
        /* Do nothing in mono */
        tmp_ch_process_method = SRLA_CH_PROCESS_METHOD_NONE;
        tmp_output_bits = code_length[0];
    } else if (header->num_channels >= 2) {
        SRLAChannelProcessMethod argmin;
        uint32_t len[4], min;

        /* Select the shortest code length among LR, MS, LS, and SR */
        SRLA_STATIC_ASSERT((SRLA_CH_PROCESS_METHOD_NONE == 0) && (SRLA_CH_PROCESS_METHOD_MS == 1)
            && (SRLA_CH_PROCESS_METHOD_LS == 2) && (SRLA_CH_PROCESS_METHOD_SR == 3));
        len[SRLA_CH_PROCESS_METHOD_NONE] = code_length[0] + code_length[1];
        len[SRLA_CH_PROCESS_METHOD_MS] = ms_code_length[0] + ms_code_length[1];
        len[SRLA_CH_PROCESS_METHOD_LS] = code_length[0] + ms_code_length[1];
        len[SRLA_CH_PROCESS_METHOD_SR] = code_length[1] + ms_code_length[1];
        min = len[SRLA_CH_PROCESS_METHOD_NONE]; argmin = SRLA_CH_PROCESS_METHOD_NONE;
        for (ch = 1; ch < 4; ch++) {
            if (min > len[ch]) {
                min = len[ch];
                argmin = (SRLAChannelProcessMethod)ch;
            }
        }

        /* Record the result */
        tmp_ch_process_method = argmin;
        tmp_output_bits = min;

        /* Replace the LR result depending on the judgment result */
        if (tmp_ch_process_method == SRLA_CH_PROCESS_METHOD_MS) {
            int32_t *tmpp;
            for (ch = 0; ch < 2; ch++) {
                memcpy(encoder->pre_emphasis[ch], encoder->ms_pre_emphasis[ch],
                    sizeof(struct SRLAPreemphasisFilter) * SRLA_NUM_PREEMPHASIS_FILTERS);
                encoder->coef_order[ch] = encoder->ms_coef_order[ch];
                encoder->coef_rshifts[ch] = encoder->ms_coef_rshifts[ch];
                memcpy(encoder->params_int[ch], encoder->ms_params_int[ch],
                    sizeof(int32_t) * encoder->ms_coef_order[ch]);
                encoder->use_sum_coef[ch] = encoder->ms_use_sum_coef[ch];
                tmpp = encoder->residual[ch];
                encoder->residual[ch] = encoder->ms_residual[ch];
                encoder->ms_residual[ch] = tmpp;
            }
        } else if ((tmp_ch_process_method == SRLA_CH_PROCESS_METHOD_LS) || (tmp_ch_process_method == SRLA_CH_PROCESS_METHOD_SR)) {
            int32_t *tmpp;
            const uint32_t src_ch = 1; /* S */
            const uint32_t dst_ch = (tmp_ch_process_method == SRLA_CH_PROCESS_METHOD_LS) ? 1 : 0;
            memcpy(encoder->pre_emphasis[dst_ch], encoder->ms_pre_emphasis[src_ch],
                sizeof(struct SRLAPreemphasisFilter) * SRLA_NUM_PREEMPHASIS_FILTERS);
            encoder->coef_order[dst_ch] = encoder->ms_coef_order[src_ch];
            encoder->coef_rshifts[dst_ch] = encoder->ms_coef_rshifts[src_ch];
            memcpy(encoder->params_int[dst_ch], encoder->ms_params_int[src_ch],
                sizeof(int32_t) * encoder->ms_coef_order[src_ch]);
            encoder->use_sum_coef[dst_ch] = encoder->ms_use_sum_coef[src_ch];
            tmpp = encoder->residual[dst_ch];
            encoder->residual[dst_ch] = encoder->ms_residual[src_ch];
            encoder->ms_residual[src_ch] = tmpp;
        }
    }

    /* Add size of multichannel processing method */
    tmp_output_bits += 2;

    /* Round up to a byte boundary */
    tmp_output_bits = SRLAUTILITY_ROUNDUP(tmp_output_bits, 8);

    /* Result output */
    (*ch_process_method) = tmp_ch_process_method;
    (*output_bits) = tmp_output_bits;

    return SRLA_APIRESULT_OK;
}

/* Compressed data block encoding */
static SRLAApiResult SRLAEncoder_EncodeCompressData(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size)
{
    uint32_t ch, tmp_code_length;
    struct BitStream writer;
    const struct SRLAHeader *header;
    SRLAChannelProcessMethod ch_process_method = SRLA_CH_PROCESS_METHOD_INVALID;

    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(encoder != NULL);
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(num_samples > 0);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(data_size > 0);
    SRLA_ASSERT(output_size != NULL);

    /* Get header */
    header = &(encoder->header);

    /* Coefficient calculation */
    if (SRLAEncoder_ComputeCoefficients(encoder,
        input, num_samples, &ch_process_method, &tmp_code_length) != SRLA_APIRESULT_OK) {
        return SRLA_APIRESULT_NG;
    }

    /* Create bit writer */
    BitWriter_Open(&writer, data, data_size);

    /* Write multichannel processing method */
    SRLA_ASSERT(ch_process_method != SRLA_CH_PROCESS_METHOD_INVALID);
    SRLA_ASSERT(ch_process_method < 4);
    BitWriter_PutBits(&writer, ch_process_method, 2);

    /* Parameter encoding */
    /* Pre-emphasis */
    for (ch = 0; ch < header->num_channels; ch++) {
        uint32_t p, uval;
        /* Pre-emphasis filter buffer */
        uval = SRLAUTILITY_SINT32_TO_UINT32(encoder->pre_emphasis[ch][0].prev);
        SRLA_ASSERT(uval < (1U << (header->bits_per_sample + 1)));
        BitWriter_PutBits(&writer, uval, header->bits_per_sample + 1);
        for (p = 0; p < SRLA_NUM_PREEMPHASIS_FILTERS; p++) {
            uval = SRLAUTILITY_SINT32_TO_UINT32(encoder->pre_emphasis[ch][p].coef);
            SRLA_ASSERT(uval < (1U << (SRLA_PREEMPHASIS_COEF_SHIFT + 1)));
            BitWriter_PutBits(&writer, uval, SRLA_PREEMPHASIS_COEF_SHIFT + 1);
        }
    }
    /* LPC coefficient order/LPC coefficient right shift amount/LPC coefficient */
    for (ch = 0; ch < header->num_channels; ch++) {
        uint32_t i, uval;
        /* LPC coefficient order */
        SRLA_ASSERT(encoder->coef_order[ch] > 0);
        SRLA_ASSERT(encoder->coef_order[ch] <= (1U << SRLA_LPC_COEFFICIENT_ORDER_BITWIDTH));
        BitWriter_PutBits(&writer, encoder->coef_order[ch] - 1, SRLA_LPC_COEFFICIENT_ORDER_BITWIDTH);
        /* LPC coefficient right shift amount */
        SRLA_ASSERT(encoder->coef_rshifts[ch] < (1U << SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH));
        BitWriter_PutBits(&writer, encoder->coef_rshifts[ch], SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH);
        /* LPC coefficients */
        BitWriter_PutBits(&writer, encoder->use_sum_coef[ch], 1); /* Flag indicating whether sum was taken and recorded */
        if (!encoder->use_sum_coef[ch]) {
            for (i = 0; i < encoder->coef_order[ch]; i++) {
                uval = SRLAUTILITY_SINT32_TO_UINT32(encoder->params_int[ch][i]);
                SRLA_ASSERT(uval < (1U << SRLA_LPC_COEFFICIENT_BITWIDTH));
                StaticHuffman_PutCode(&encoder->param_codes, &writer, uval);
            }
        } else {
            uval = SRLAUTILITY_SINT32_TO_UINT32(encoder->params_int[ch][0]);
            SRLA_ASSERT(uval < (1U << SRLA_LPC_COEFFICIENT_BITWIDTH));
            StaticHuffman_PutCode(&encoder->param_codes, &writer, uval);
            for (i = 1; i < encoder->coef_order[ch]; i++) {
                const int32_t summed = encoder->params_int[ch][i] + encoder->params_int[ch][i - 1];
                uval = SRLAUTILITY_SINT32_TO_UINT32(summed);
                SRLA_ASSERT(uval < (1U << SRLA_LPC_COEFFICIENT_BITWIDTH));
                StaticHuffman_PutCode(&encoder->sum_param_codes, &writer, uval);
            }
        }
    }

    /* Residual coding */
    for (ch = 0; ch < header->num_channels; ch++) {
        SRLACoder_Encode(encoder->coder, &writer, encoder->residual[ch], num_samples);
    }

    /* Align to byte boundary */
    BitStream_Flush(&writer);

    /* Get write size */
    BitStream_Tell(&writer, (int32_t *)output_size);

    /* Discard bit writer */
    BitStream_Close(&writer);

    return SRLA_APIRESULT_OK;
}

/* Silence data block encoding */
static SRLAApiResult SRLAEncoder_EncodeSilentData(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size)
{
    /* Since this is an internal function, invalid arguments are dropped with an assertion */
    SRLA_ASSERT(encoder != NULL);
    SRLA_ASSERT(input != NULL);
    SRLA_ASSERT(num_samples > 0);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(data_size > 0);
    SRLA_ASSERT(output_size != NULL);

    /* No data size */
    (*output_size) = 0;
    return SRLA_APIRESULT_OK;
}

/* Calculate the size of a single data block */
SRLAApiResult SRLAEncoder_ComputeBlockSize(
    struct SRLAEncoder *encoder, const int32_t *const *input, uint32_t num_samples,
    uint32_t *output_size)
{
    uint8_t *data_ptr;
    const struct SRLAHeader *header;
    SRLABlockDataType block_type;
    uint32_t tmp_block_size;

    /* Argument check */
    if ((encoder == NULL) || (input == NULL) || (num_samples == 0)
        || (output_size == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }
    header = &(encoder->header);

    /* No parameters set */
    if (encoder->set_parameter != 1) {
        return SRLA_APIRESULT_PARAMETER_NOT_SET;
    }

    /* Check the number of encoded samples */
    if (num_samples > header->max_num_samples_per_block) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Determine compression method */
    block_type = SRLAEncoder_DecideBlockDataType(encoder, input, num_samples);
    SRLA_ASSERT(block_type != SRLA_BLOCK_DATA_TYPE_INVALID);

COMPUTE_BLOCK_SIZE_START:
    /* Block header size */
    tmp_block_size = 11;

    switch (block_type) {
    case SRLA_BLOCK_DATA_TYPE_RAWDATA:
        tmp_block_size += (header->bits_per_sample * num_samples * header->num_channels) / 8;
        break;
    case SRLA_BLOCK_DATA_TYPE_COMPRESSDATA:
    {
        SRLAApiResult ret;
        uint32_t compress_data_size;
        SRLAChannelProcessMethod dummy;
        /* Code length calculation */
        if ((ret = SRLAEncoder_ComputeCoefficients(
            encoder, input, num_samples, &dummy, &compress_data_size)) != SRLA_APIRESULT_OK) {
            return ret;
        }
        SRLA_ASSERT(compress_data_size % 8 == 0);
        /* Switch to raw data block when data increases as a result of encoding */
        if (compress_data_size >= (header->bits_per_sample * num_samples * header->num_channels)) {
            block_type = SRLA_BLOCK_DATA_TYPE_RAWDATA;
            goto COMPUTE_BLOCK_SIZE_START;
        }
        /* Convert to bytes and add */
        tmp_block_size += compress_data_size / 8;
    }
        break;
    case SRLA_BLOCK_DATA_TYPE_SILENT:
        /* 0 bytes */
        break;
    default:
        SRLA_ASSERT(0);
    }

    /* Result output */
    (*output_size) = tmp_block_size;

    return SRLA_APIRESULT_OK;
}

/* Single data block encoding */
SRLAApiResult SRLAEncoder_EncodeBlock(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size)
{
    uint8_t *data_ptr;
    const struct SRLAHeader *header;
    SRLABlockDataType block_type;
    SRLAApiResult ret;
    uint32_t block_header_size, block_data_size;

    /* Argument check */
    if ((encoder == NULL) || (input == NULL) || (num_samples == 0)
            || (data == NULL) || (data_size == 0) || (output_size == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }
    header = &(encoder->header);

    /* No parameters set */
    if (encoder->set_parameter != 1) {
        return SRLA_APIRESULT_PARAMETER_NOT_SET;
    }

    /* Check the number of encoded samples */
    if (num_samples > header->max_num_samples_per_block) {
        return SRLA_APIRESULT_INSUFFICIENT_BUFFER;
    }

    /* Determine compression method */
    block_type = SRLAEncoder_DecideBlockDataType(encoder, input, num_samples);
    SRLA_ASSERT(block_type != SRLA_BLOCK_DATA_TYPE_INVALID);

ENCODING_BLOCK_START:
    /* Encode the block header */
    data_ptr = data;
    /* Synchronization code at the beginning of a block */
    ByteArray_PutUint16BE(data_ptr, SRLA_BLOCK_SYNC_CODE);
    /* Block size: fill with temporary value */
    ByteArray_PutUint32BE(data_ptr, 0);
    /* Block checksum: fill with a temporary value */
    ByteArray_PutUint16BE(data_ptr, 0);
    /* Block data type */
    ByteArray_PutUint8(data_ptr, block_type);
    /* Number of samples per block channel */
    ByteArray_PutUint16BE(data_ptr, num_samples);
    /* Block header size */
    block_header_size = (uint32_t)(data_ptr - data);

    /* Data part encoding */
    /* Call different encoding functions depending on the method */
    switch (block_type) {
    case SRLA_BLOCK_DATA_TYPE_RAWDATA:
        ret = SRLAEncoder_EncodeRawData(encoder, input, num_samples,
                data_ptr, data_size - block_header_size, &block_data_size);
        break;
    case SRLA_BLOCK_DATA_TYPE_COMPRESSDATA:
        ret = SRLAEncoder_EncodeCompressData(encoder, input, num_samples,
                data_ptr, data_size - block_header_size, &block_data_size);
        /* Switch to raw data block when data increases as a result of encoding */
        if ((8 * block_data_size) >= (header->bits_per_sample * num_samples * header->num_channels)) {
            block_type = SRLA_BLOCK_DATA_TYPE_RAWDATA;
            goto ENCODING_BLOCK_START;
        }
        break;
    case SRLA_BLOCK_DATA_TYPE_SILENT:
        ret = SRLAEncoder_EncodeSilentData(encoder, input, num_samples,
                data_ptr, data_size - block_header_size, &block_data_size);
        break;
    default:
        ret = SRLA_APIRESULT_INVALID_FORMAT;
        break;
    }

    /* Encoding failed */
    if (ret != SRLA_APIRESULT_OK) {
        return ret;
    }

    /*
/* Write block size:
* Checksum (2 bytes) + number of samples per block channel (2 bytes) + block data type (1 byte) */
*/
    ByteArray_WriteUint32BE(&data[2], block_data_size + 5);

    /* Calculate and write the checksum after the checksum area */
    {
        /* Add the number of samples per block channel (2 bytes) + block data type (1 byte) */
        const uint16_t checksum = SRLAUtility_CalculateFletcher16CheckSum(&data[8], block_data_size + 3);
        ByteArray_WriteUint16BE(&data[6], checksum);
    }

    /* Output size */
    (*output_size) = block_header_size + block_data_size;

    /* Encoding successful */
    return SRLA_APIRESULT_OK;
}

/* Encoding including optimal block division search */
SRLAApiResult SRLAEncoder_EncodeOptimalPartitionedBlock(
    struct SRLAEncoder *encoder,
    const int32_t *const *input, uint32_t num_samples,
    uint8_t *data, uint32_t data_size, uint32_t *output_size)
{
    SRLAApiResult ret;
    uint32_t num_partitions, part, ch;
    uint32_t progress, write_offset, tmp_output_size;

    /* Argument check */
    if ((encoder == NULL) || (input == NULL)
        || (data == NULL) || (output_size == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* No parameters set */
    if (encoder->set_parameter != 1) {
        return SRLA_APIRESULT_PARAMETER_NOT_SET;
    }

    /* Finding the best block division */
    if (SRLAEncoder_SearchOptimalBlockPartitions(
        encoder, input, num_samples, encoder->min_num_samples_per_block,
        &num_partitions, encoder->partitions_buffer) != SRLA_ERROR_OK) {
        return SRLA_APIRESULT_NG;
    }

    /* Encode according to division */
    progress = write_offset = 0;
    for (part = 0; part < num_partitions; part++) {
        const uint32_t num_block_samples = encoder->partitions_buffer[part];
        const int32_t *input_ptr[SRLA_MAX_NUM_CHANNELS];
        for (ch = 0; ch < encoder->header.num_channels; ch++) {
            input_ptr[ch] = &input[ch][progress];
        }
        if ((ret = SRLAEncoder_EncodeBlock(encoder,
                input_ptr, num_block_samples, data + write_offset, data_size - write_offset,
                &tmp_output_size)) != SRLA_APIRESULT_OK) {
            return ret;
        }
        write_offset += tmp_output_size;
        progress += num_block_samples;
        SRLA_ASSERT(write_offset <= data_size);
        SRLA_ASSERT(progress <= num_samples);
    }
    SRLA_ASSERT(progress == num_samples);

    /* Successful completion */
    (*output_size) = write_offset;
    return SRLA_APIRESULT_OK;
}

/* Encode the entire file, including the header */
SRLAApiResult SRLAEncoder_EncodeWhole(
        struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size, uint8_t variable_block)
{
    SRLAApiResult ret;
    uint32_t progress, ch, write_size, write_offset, num_encode_samples;
    uint8_t *data_pos;
    const int32_t *input_ptr[SRLA_MAX_NUM_CHANNELS];
    const struct SRLAHeader *header;
    SRLAApiResult (*encode_function)(struct SRLAEncoder *encoder,
        const int32_t *const *input, uint32_t num_samples,
        uint8_t *data, uint32_t data_size, uint32_t *output_size);

    /* Argument check */
    if ((encoder == NULL) || (input == NULL)
            || (data == NULL) || (output_size == NULL)) {
        return SRLA_APIRESULT_INVALID_ARGUMENT;
    }

    /* No parameters set */
    if (encoder->set_parameter != 1) {
        return SRLA_APIRESULT_PARAMETER_NOT_SET;
    }

    /* Select the encoding function */
    encode_function = (variable_block == 1) ? SRLAEncoder_EncodeOptimalPartitionedBlock : SRLAEncoder_EncodeBlock;

    /* Get the write position */
    data_pos = data;

    /* Header encoding */
    encoder->header.num_samples = num_samples;
    if ((ret = SRLAEncoder_EncodeHeader(&(encoder->header), data_pos, data_size))
            != SRLA_APIRESULT_OK) {
        return ret;
    }
    header = &(encoder->header);

    /* Initialize progress status */
    progress = 0;
    write_offset = SRLA_HEADER_SIZE;
    data_pos = data + SRLA_HEADER_SIZE;

    /* Encode blocks in chronological order */
    while (progress < num_samples) {
        /* Determine the number of samples to encode */
        num_encode_samples
            = SRLAUTILITY_MIN(header->max_num_samples_per_block, num_samples - progress);

        /* Set sample reference position */
        for (ch = 0; ch < header->num_channels; ch++) {
            input_ptr[ch] = &input[ch][progress];
        }

        /* Block encoding */
        if ((ret = encode_function(encoder,
                input_ptr, num_encode_samples,
                data_pos, data_size - write_offset, &write_size)) != SRLA_APIRESULT_OK) {
            return ret;
        }

        /* Progress update */
        data_pos += write_size;
        write_offset += write_size;
        progress += num_encode_samples;
        SRLA_ASSERT(write_offset <= data_size);
    }

    /* Successful completion */
    (*output_size) = write_offset;
    return SRLA_APIRESULT_OK;
}
