#ifndef LPC_H_INCLUDED
#define LPC_H_INCLUDED

#include <stdint.h>

/* API result type */
typedef enum LPCApiResultTag {
    LPC_APIRESULT_OK = 0,                 /* OK */
    LPC_APIRESULT_NG,                     /* Unclassifiable error */
    LPC_APIRESULT_INVALID_ARGUMENT,       /* Invalid argument */
    LPC_APIRESULT_EXCEED_MAX_ORDER,       /* Maximum degree exceeded */
    LPC_APIRESULT_EXCEED_MAX_NUM_SAMPLES, /* Maximum input samples exceeded */
    LPC_APIRESULT_FAILED_TO_CALCULATION   /* Calculation failed */
} LPCApiResult;

/* Window function type */
typedef enum LPCWindowTypeTag {
    LPC_WINDOWTYPE_RECTANGULAR = 0, /* Rectangular window */
    LPC_WINDOWTYPE_SIN,             /* Sign window */
    LPC_WINDOWTYPE_WELCH            /* Welch window */
} LPCWindowType;

/* LPC coefficient calculation handle */
struct LPCCalculator;

/* Initialization config */
struct LPCCalculatorConfig {
    uint32_t max_order;        /* Maximum degree */
    uint32_t max_num_samples;  /* Maximum number of input samples */
};

#ifdef __cplusplus
extern "C" {
#endif

/* Calculate the work size of the LPC coefficient calculation handle */
int32_t LPCCalculator_CalculateWorkSize(const struct LPCCalculatorConfig *config);

/* Create LPC coefficient calculation handle */
struct LPCCalculator *LPCCalculator_Create(const struct LPCCalculatorConfig *config, void *work, int32_t work_size);

/* Discard LPC coefficient calculation handle */
void LPCCalculator_Destroy(struct LPCCalculator *lpcc);

/* Calculate LPC coefficients using Levinson-Durbin recursive calculation */
LPCApiResult LPCCalculator_CalculateLPCCoefficients(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    LPCWindowType window_type, double regular_term);

/* Calculate all LPC coefficients up to a given order using Levinson-Durbin recursion (double precision) */
/* error_vars is calculated from the 0th order error variance (variance) to the max_coef_order order variance, so the size of error_vars needs to be max_coef_order+1 */
LPCApiResult LPCCalculator_CalculateMultipleLPCCoefficients(
    struct LPCCalculator* lpcc,
    const double* data, uint32_t num_samples, double **lpc_coefs, double *error_vars, uint32_t max_coef_order,
    LPCWindowType window_type, double regular_term);

/* Calculate LPC coefficients using auxiliary function method (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsAF(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    uint32_t max_num_iteration, LPCWindowType window_type, double regular_term);

/* Calculate LPC coefficients using Burg method (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsBurg(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order);

/* Calculate LPC coefficients from SVR (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsSVR(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    uint32_t max_num_iteration, LPCWindowType window_type,
    double regular_term, const double *margin_list, uint32_t margin_list_size);

/* Calculate the estimated code length per sample from the input data */
LPCApiResult LPCCalculator_EstimateCodeLength(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, uint32_t bits_per_sample,
    uint32_t coef_order, double *length_per_sample_bits, LPCWindowType window_type);

/* Calculate the MDL (Minimum Description Length) */
LPCApiResult LPCCalculator_CalculateMDL(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, uint32_t coef_order, double *mdl,
    LPCWindowType window_type);

/* Convert LPC coefficients to PARCOR coefficients and quantize them */
LPCApiResult LPC_QuantizeCoefficientsAsPARCOR(
    struct LPCCalculator *lpcc,
    const double *lpc_coef, uint32_t coef_order, uint32_t nbits_precision, int32_t *int_coef);

/* Integer quantization of LPC coefficients */
LPCApiResult LPC_QuantizeCoefficients(
    const double *double_coef, uint32_t coef_order, uint32_t nbits_precision, uint32_t max_bits,
    int32_t *int_coef, uint32_t *coef_rshift);

/* Prediction/error output by LPC coefficients */
LPCApiResult LPC_Predict(
    const int32_t *data, uint32_t num_samples,
    const int32_t *coef, uint32_t coef_order, int32_t *residual, uint32_t coef_rshift);

/* Synthesis by LPC coefficients (in-place) */
LPCApiResult LPC_Synthesize(
    int32_t *data, uint32_t num_samples,
    const int32_t *coef, uint32_t coef_order, uint32_t coef_rshift);

#ifdef __cplusplus
}
#endif

#endif /* LPC_H_INCLUDED */
