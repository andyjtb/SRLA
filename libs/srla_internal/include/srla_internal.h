#ifndef SRLA_INTERNAL_H_INCLUDED
#define SRLA_INTERNAL_H_INCLUDED

#include "srla.h"
#include "srla_stdint.h"
#include "static_huffman.h"

/* Memory alignment of this library */
#define SRLA_MEMORY_ALIGNMENT 16
/* Synchronization code at the beginning of a block */
#define SRLA_BLOCK_SYNC_CODE 0xFFFF

/* Internal encoding parameters */
/* Pre-emphasis coefficient shift amount */
#define SRLA_PREEMPHASIS_COEF_SHIFT 4
/* Number of times to apply pre-emphasis filter */
#define SRLA_NUM_PREEMPHASIS_FILTERS 2
/* Bit width of LPC coefficients */
#define SRLA_LPC_COEFFICIENT_BITWIDTH 8
/* Bit width of LPC coefficient right shift amount */
#define SRLA_RSHIFT_LPC_COEFFICIENT_BITWIDTH 4
/* (LPC coefficient order - 1) bit width */
#define SRLA_LPC_COEFFICIENT_ORDER_BITWIDTH 8
/* Ridge regularization parameter for LPC */
#define SRLA_LPC_RIDGE_REGULARIZATION_PARAMETER 1e-5

/* Assert macro */
#ifdef NDEBUG
/* Explicitly avoid unused variable warnings */
#define SRLA_ASSERT(condition) ((void)(condition))
#else
#include <assert.h>
#define SRLA_ASSERT(condition) assert(condition)
#endif

/* static assert macro */
#define SRLA_STATIC_ASSERT(expr) extern void assertion_failed(char dummy[(expr) ? 1 : -1])

/* Block data type */
typedef enum SRLABlockDataTypeTag {
    SRLA_BLOCK_DATA_TYPE_COMPRESSDATA  = 0, /* Compressed data */
    SRLA_BLOCK_DATA_TYPE_SILENT        = 1, /* Silence data */
    SRLA_BLOCK_DATA_TYPE_RAWDATA       = 2, /* Raw data */
    SRLA_BLOCK_DATA_TYPE_INVALID       = 3  /* invalid           */
} SRLABlockDataType;

/* How to determine multi-channel processing */
typedef enum SRLAChannelProcessMethodTacticsTag {
    SRLA_CH_PROCESS_METHOD_TACTICS_NONE = 0, /* Do nothing */
    SRLA_CH_PROCESS_METHOD_TACTICS_MS_FIXED, /* Always select stereo MS processing */
    SRLA_CH_PROCESS_METHOD_TACTICS_ADAPTIVE, /* Adaptively select LR, LS, RS, MS */
    SRLA_CH_PROCESS_METHOD_TACTICS_INVALID   /* Invalid value */
} SRLAChannelProcessMethodTactics;

/* Multi-channel processing method */
typedef enum SRLAChannelProcessMethodTag {
    SRLA_CH_PROCESS_METHOD_NONE = 0, /* Do nothing */
    SRLA_CH_PROCESS_METHOD_MS = 1, /* Stereo MS processing */
    SRLA_CH_PROCESS_METHOD_LS = 2, /* Stereo LS processing */
    SRLA_CH_PROCESS_METHOD_SR = 3, /* Stereo SR processing */
    SRLA_CH_PROCESS_METHOD_INVALID /* Invalid value */
} SRLAChannelProcessMethod;

/* How to determine the LPC order */
typedef enum SRLAChannelLPCOrderDecisionTacticsTag {
    SRLA_LPC_ORDER_DECISION_TACTICS_MAX_FIXED = 0, /* Always choose the highest order */
    SRLA_LPC_ORDER_DECISION_TACTICS_BRUTEFORCE_SEARCH, /* Simple exhaustive search */
    SRLA_LPC_ORDER_DECISION_TACTICS_BRUTEFORCE_ESTIMATION, /* Exhaustive search by estimating residual variance */
    SRLA_LPC_ORDER_DECISION_TACTICS_INVALID  /* Invalid value */
} SRLAChannelLPCOrderDecisionTactics;

/* Internal error type */
typedef enum SRLAErrorTag {
    SRLA_ERROR_OK = 0, /* OK */
    SRLA_ERROR_NG, /* Unclassifiable failure */
    SRLA_ERROR_INVALID_ARGUMENT, /* Invalid argument */
    SRLA_ERROR_INVALID_FORMAT, /* Invalid format */
    SRLA_ERROR_INSUFFICIENT_BUFFER, /* Buffer size is insufficient */
    SRLA_ERROR_INSUFFICIENT_DATA /* Data size is insufficient */
} SRLAError;

/* Parameter presets */
struct SRLAParameterPreset {
    uint32_t max_num_parameters; /* Maximum number of parameters */
    SRLAChannelProcessMethodTactics ch_process_method_tactics; /* How to determine multi-channel processing */
    SRLAChannelLPCOrderDecisionTactics lpc_order_tactics; /* LPC order determination method */
    uint32_t svr_max_num_iterations; /* Maximum number of SVR iterations */
    const double *margin_list; /* margin list */
    uint32_t margin_list_size; /* margin list size */
};

#ifdef __cplusplus
extern "C" {
#endif

/* Parameter preset array */
extern const struct SRLAParameterPreset g_srla_parameter_preset[];

/* Get the Huffman tree for the parameter codes */
const struct StaticHuffmanTree* SRLA_GetParameterHuffmanTree(void);

/* Get the Huffman tree for the summed parameter codes */
const struct StaticHuffmanTree *SRLA_GetSumParameterHuffmanTree(void);


#ifdef __cplusplus
}
#endif

#endif /* SRLA_INTERNAL_H_INCLUDED */
