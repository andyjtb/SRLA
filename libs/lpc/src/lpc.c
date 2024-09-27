#include "lpc.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

#include "fft.h"

/* memory alignment */
#define LPC_ALIGNMENT 16
/* Pi */
#define LPC_PI 3.1415926535897932384626433832795029
/* Minimum absolute residual value */
#define LPCAF_RESIDUAL_EPSILON 1e-6

/* Round up to the next multiple of n */
#define LPC_ROUNDUP(val, n) ((((val) + ((n) - 1)) / (n)) * (n))
/* Sign function */
#define LPC_SIGN(val) (((val) > 0) - ((val) < 0))
/* Get the maximum value */
#define LPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
/* Get absolute value */
#define LPC_ABS(val) (((val) > 0) ? (val) : -(val))
/* Soft threshold operator */
#define LPC_SOFT_THRESHOLD(in, epsilon) (LPC_SIGN(in) * LPC_MAX(LPC_ABS(in) - (epsilon), 0.0))

/* Internal error type */
typedef enum LPCErrorTag {
    LPC_ERROR_OK = 0,
    LPC_ERROR_NG,
    LPC_ERROR_SINGULAR_MATRIX,
    LPC_ERROR_INVALID_ARGUMENT
} LPCError;

/* LPC calculation handle */
struct LPCCalculator {
    uint32_t max_order; /* Maximum degree */
    uint32_t max_num_buffer_samples; /* Maximum number of buffer samples */
    /* All internal calculation results are stored as doubles to ensure accuracy */
    /* If you increase the number of samples with float, the output calculation result will be nan due to an error in the sample autocorrelation value. */
    double **a_vecs; /* Coefficient vector for each order */
    double *u_vec; /* Calculation vector 3 */
    double *v_vec; /* Calculation vector 4 */
    double **r_mat; /* Matrix used in auxiliary function method/Burg method ((max_order + 1) order) */
    double *auto_corr; /* Sample autocorrelation */
    double *parcor_coef; /* PARCOR coefficient vector */
    double *error_vars; /* Residual variance */
    double *buffer; /* Input signal buffer area */
    double *work_buffer; /* Calculation buffer */
    uint8_t alloced_by_own; /* Did you allocate the space yourself? */
    void *work; /* Work area start pointer */
};

/* round function (not defined in C89) */
static double LPC_Round(double d)
{
    return (d >= 0.0) ? floor(d + 0.5) : -floor(-d + 0.5);
}

/* log2 function (not defined in C89) */
static double LPC_Log2(double d)
{
#define INV_LOGE2 (1.4426950408889634)  /* 1 / log(2) */
    return log(d) * INV_LOGE2;
#undef INV_LOGE2
}

/* Round up to a power of 2 */
static uint32_t LPC_RoundUp2Powered(uint32_t val)
{
    /* Hacker's Fun Reference */
    val--;
    val |= val >> 1;
    val |= val >> 2;
    val |= val >> 4;
    val |= val >> 8;
    val |= val >> 16;
    return val + 1;
}

/* Calculate the work size of the LPC coefficient calculation handle */
int32_t LPCCalculator_CalculateWorkSize(const struct LPCCalculatorConfig *config)
{
    int32_t work_size;

    /* Argument check */
    if (config == NULL) {
        return -1;
    }

    work_size = sizeof(struct LPCCalculator) + LPC_ALIGNMENT;
    /* Area to be used for a_vec */
    work_size += (int32_t)(sizeof(double *) * (config->max_order + 1));
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1) * (config->max_order + 2));
    /* Area for u, v vectors */
    work_size += (int32_t)(sizeof(double) * (config->max_order + 2) * 2);
    /* domain of sample autocorrelation */
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1));
    /* Domain of LPC coefficient vector */
    work_size += (int32_t)(sizeof(double *) * (config->max_order + 1));
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1) * (config->max_order + 1));
    /* Domain of PARCOR coefficient vector */
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1));
    /* domain of residual variance */
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1));
    /* Matrix domain used for auxiliary function method */
    work_size += (int32_t)(sizeof(double *) * (config->max_order + 1));
    work_size += (int32_t)(sizeof(double) * (config->max_order + 1) * (config->max_order + 1));
    /* Input signal buffer area */
    work_size += (int32_t)(sizeof(double) * LPC_RoundUp2Powered(config->max_num_samples));
    /* Buffer area for calculation */
    work_size += (int32_t)(sizeof(double) * LPC_RoundUp2Powered(config->max_num_samples));

    return work_size;
}

/* Create LPC coefficient calculation handle */
struct LPCCalculator* LPCCalculator_Create(const struct LPCCalculatorConfig *config, void *work, int32_t work_size)
{
    struct LPCCalculator *lpcc;
    uint8_t *work_ptr;
    uint8_t tmp_alloc_by_own = 0;

    /* Allocate work area yourself */
    if ((work == NULL) && (work_size == 0)) {
        if ((work_size = LPCCalculator_CalculateWorkSize(config)) < 0) {
            return NULL;
        }
        work = malloc((uint32_t)work_size);
        tmp_alloc_by_own = 1;
    }

    /* Argument check */
    if ((config == NULL) || (work == NULL)
            || (work_size < LPCCalculator_CalculateWorkSize(config))
            || (config->max_order == 0) || (config->max_num_samples == 0)) {
        if (tmp_alloc_by_own == 1) {
            free(work);
        }
        return NULL;
    }

    /* Get work area */
    work_ptr = (uint8_t *)work;

    /* Allocate handle area */
    work_ptr = (uint8_t *)LPC_ROUNDUP((uintptr_t)work_ptr, LPC_ALIGNMENT);
    lpcc = (struct LPCCalculator *)work_ptr;
    work_ptr += sizeof(struct LPCCalculator);

    /* Set handle members */
    lpcc->max_order = config->max_order;
    lpcc->max_num_buffer_samples = config->max_num_samples;
    lpcc->work = work;
    lpcc->alloced_by_own = tmp_alloc_by_own;

    /* Allocation of vector space for calculation */
    {
        uint32_t ord;
        lpcc->a_vecs = (double **)work_ptr;
        work_ptr += sizeof(double *) * (config->max_order + 1);
        for (ord = 0; ord < config->max_order + 1; ord++) {
            lpcc->a_vecs[ord] = (double *)work_ptr;
            work_ptr += sizeof(double) * (config->max_order + 2); /* max_order+2 including a_0, a_k+1 */
        }
    }
    lpcc->u_vec = (double *)work_ptr;
    work_ptr += sizeof(double) * (config->max_order + 2);
    lpcc->v_vec = (double *)work_ptr;
    work_ptr += sizeof(double) * (config->max_order + 2);

    /* Sample autocorrelation domain allocation */
    lpcc->auto_corr = (double *)work_ptr;
    work_ptr += sizeof(double) * (config->max_order + 1);

    /* Space allocation for PARCOR coefficient vector */
    lpcc->parcor_coef = (double *)work_ptr;
    work_ptr += sizeof(double) * (config->max_order + 1);

    /* Residual variance domain allocation */
    lpcc->error_vars = (double *)work_ptr;
    work_ptr += sizeof(double) * (config->max_order + 1);

    /* Matrix domain used in auxiliary function method/Burg method */
    {
        uint32_t ord;
        lpcc->r_mat = (double **)work_ptr;
        work_ptr += sizeof(double *) * (config->max_order + 1);
        for (ord = 0; ord < config->max_order + 1; ord++) {
            lpcc->r_mat[ord] = (double *)work_ptr;
            work_ptr += sizeof(double) * (config->max_order + 1);
        }
    }

    /* Input signal buffer area */
    lpcc->buffer = (double *)work_ptr;
    work_ptr += sizeof(double) * LPC_RoundUp2Powered(config->max_num_samples);

    /* Calculation buffer area */
    lpcc->work_buffer = (double *)work_ptr;
    work_ptr += sizeof(double) * LPC_RoundUp2Powered(config->max_num_samples);

    /* Buffer overflow check */
    assert((work_ptr - (uint8_t *)work) <= work_size);

    return lpcc;
}

/* Discard LPC coefficient calculation handle */
void LPCCalculator_Destroy(struct LPCCalculator *lpcc)
{
    if (lpcc != NULL) {
        /* Release the work area if it was previously allocated */
        if (lpcc->alloced_by_own == 1) {
            free(lpcc->work);
        }
    }
}

/* Apply window function */
static LPCError LPC_ApplyWindow(
    LPCWindowType window_type, const double *input, uint32_t num_samples, double *output)
{
    /* Argument check */
    if (input == NULL || output == NULL) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    switch (window_type) {
    case LPC_WINDOWTYPE_RECTANGULAR:
        memcpy(output, input, sizeof(double) * num_samples);
        break;
    case LPC_WINDOWTYPE_SIN:
        {
            uint32_t smpl;
            for (smpl = 0; smpl < num_samples; smpl++) {
                output[smpl] = input[smpl] * sin((LPC_PI * smpl) / (num_samples - 1));
            }
        }
        break;
    case LPC_WINDOWTYPE_WELCH:
        {
            uint32_t smpl;
            const double divisor = 4.0 * pow(num_samples - 1, -2.0);
            for (smpl = 0; smpl < (num_samples >> 1); smpl++) {
                const double weight = divisor * smpl * (num_samples - 1 - smpl);
                output[smpl] = input[smpl] * weight;
                output[num_samples - smpl - 1] = input[num_samples - smpl - 1] * weight;
            }
        }
        break;
    default:
        return LPC_ERROR_NG;
    }

    return LPC_ERROR_OK;
}

/* Calculate the inverse of the sum of squares of the window */
static double LPC_ComputeWindowInverseSquaredSum(LPCWindowType window_type, const uint32_t num_samples)
{
    switch (window_type) {
    case LPC_WINDOWTYPE_RECTANGULAR:
        return 1.0;
    case LPC_WINDOWTYPE_WELCH:
    {
        const double n = num_samples - 1;
        return (15 * (n - 1) * (n - 1) * (n - 1)) / (8 * n * (n - 2) * (n * n - 2 * n + 2));
    }
    default:
        assert(0);
    }

    return 1.0;
}

/* (sample) autocorrelation calculation */
static LPCError LPC_CalculateAutoCorrelation(
    const double *data, uint32_t num_samples, double *auto_corr, uint32_t order)
{
    uint32_t i, lag;
    double tmp;

    assert(num_samples >= order);

    /* Argument check */
    if (data == NULL || auto_corr == NULL) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /* Coefficient initialization */
    for (lag = 0; lag < order; lag++) {
        auto_corr[lag] = 0.0;
    }

    /* Calculate the autocorrelation coefficient by focusing on the lag of the data instead of the order */
    for (i = 0; i <= num_samples - order; i++) {
        tmp = data[i];
        /* Take the sum of products of data with the same lag */
        for (lag = 0; lag < order; lag++) {
            auto_corr[lag] += tmp * data[i + lag];
        }
    }
    for (; i < num_samples; i++) {
        tmp = data[i];
        for (lag = 0; lag < num_samples - i; lag++) {
            auto_corr[lag] += tmp * data[i + lag];
        }
    }

    return LPC_ERROR_OK;
}

/* Calculates (sample) autocorrelation using FFT. The contents of data_buffer are destroyed. */
static LPCError LPC_CalculateAutoCorrelationByFFT(
    double *data_buffer, double *work_buffer, uint32_t num_buffer_samples, uint32_t num_samples, double *auto_corr, uint32_t order)
{
    uint32_t i;
    uint32_t fft_size;
    const double norm_factor = 2.0 / num_samples;

    assert(num_samples >= order);
    assert(num_buffer_samples >= num_samples);

    /* Argument check */
    if ((data_buffer == NULL) || (auto_corr == NULL)) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /* Determine FFT size */
    fft_size = LPC_RoundUp2Powered(num_samples);
    assert(num_buffer_samples >= fft_size);

    /* Fill the latter half with 0 */
    for (i = num_samples; i < fft_size; i++) {
        data_buffer[i] = 0.0;
    }

    /* FFT */
    FFT_RealFFT(fft_size, -1, data_buffer, work_buffer);

    /* Calculate the square of the complex absolute value */
    data_buffer[0] *= data_buffer[0];
    data_buffer[1] *= data_buffer[1];
    for (i = 2; i < fft_size; i += 2) {
        const double real = data_buffer[i + 0];
        const double imag = data_buffer[i + 1];
        data_buffer[i + 0] = real * real + imag * imag;
        data_buffer[i + 1] = 0.0;
    }

    /* IFFT */
    FFT_RealFFT(fft_size, 1, data_buffer, work_buffer);

    /* Result set returning normalization constants */
    for (i = 0; i < order; i++) {
        auto_corr[i] = data_buffer[i] * norm_factor;
    }

    return LPC_ERROR_OK;
}

/* Levinson-Durbin recursion */
static LPCError LPC_LevinsonDurbinRecursion(struct LPCCalculator *lpcc,
    const double *auto_corr, uint32_t coef_order, double *parcor_coef, double *error_vars)
{
    uint32_t k, i;
    double gamma; /* reflection coefficient */

    /* Copy pointer to auto variable */
    double **a_vecs = lpcc->a_vecs;

    /* Argument check */
    if ((lpcc == NULL) || (auto_corr == NULL) || (parcor_coef == NULL)) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /*
/* If the zeroth-order autocorrelation (sum of squares of the signal) is small
* => Predict a silent output system with all coefficients set to 0 */
*/
    if (fabs(auto_corr[0]) < FLT_EPSILON) {
        for (i = 0; i < coef_order + 1; i++) {
            parcor_coef[i] = 0.0;
            error_vars[i] = auto_corr[0]; /* Residual variance is the same as the input */
        }
        for (k = 0; k < coef_order; k++) {
            for (i = 0; i < coef_order + 2; i++) {
                lpcc->a_vecs[k][i] = 0.0;
            }
        }
        return LPC_ERROR_OK;
    }

    /* Set the coefficients for the first step */
    a_vecs[0][0] = 1.0;
    error_vars[0] = auto_corr[0];
    a_vecs[0][1] = - auto_corr[1] / auto_corr[0];
    a_vecs[0][2] = 0.0;
    parcor_coef[0] = auto_corr[1] / error_vars[0];
    error_vars[1] = error_vars[0] + auto_corr[1] * a_vecs[0][1];

    /* Recursive processing */
    for (k = 1; k < coef_order; k++) {
        const double *a_vec = a_vecs[k - 1];

        gamma = 0.0;
        for (i = 0; i < k + 1; i++) {
            gamma += a_vec[i] * auto_corr[k + 1 - i];
        }
        gamma /= -error_vars[k];
        error_vars[k + 1] = error_vars[k] * (1.0 - gamma * gamma);
        /* Error variance (power) is non-negative */
        assert(error_vars[k + 1] >= 0.0);

        /* Update coefficients */
        for (i = 0; i < k + 2; i++) {
            a_vecs[k][i] = a_vec[i] + gamma * a_vec[k + 1 - i];
        }
        a_vecs[k][k + 2] = 0.0;
        /* PARCOR coefficients are the sign-negated reflection coefficients */
        parcor_coef[k] = -gamma;
        /* The absolute value of the PARCOR coefficient is less than 1 (convergence condition) */
        assert(fabs(gamma) < 1.0);
    }

    return LPC_ERROR_OK;
}

/* Common functions for coefficient calculation */
static LPCError LPC_CalculateCoef(
    struct LPCCalculator *lpcc, const double *data, uint32_t num_samples, uint32_t coef_order,
    LPCWindowType window_type, double regular_term)
{
    /* Argument check */
    if (lpcc == NULL) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /* Apply window function */
    if (LPC_ApplyWindow(window_type, data, num_samples, lpcc->buffer) != LPC_ERROR_OK) {
        return LPC_ERROR_NG;
    }

    /* Calculate the autocorrelation */
#if 0
    if (LPC_CalculateAutoCorrelation(
            lpcc->buffer, num_samples, lpcc->auto_corr, coef_order + 1) != LPC_ERROR_OK) {
        return LPC_ERROR_NG;
    }
#else
    if (LPC_CalculateAutoCorrelationByFFT(
            lpcc->buffer, lpcc->work_buffer, lpcc->max_num_buffer_samples,
            num_samples, lpcc->auto_corr, coef_order + 1) != LPC_ERROR_OK) {
        return LPC_ERROR_NG;
    }
#endif

    /*
/* When the number of input samples is small, the coefficients often diverge.
* => Treat it as silent data and set all coefficients to 0.*/
*/
    if (num_samples < coef_order) {
        uint32_t i;
        for (i = 0; i < coef_order + 1; i++) {
            lpcc->parcor_coef[i] = 0.0;
        }
        return LPC_ERROR_OK;
    }

    /* Emphasize zero-order correlation (Ridge regularization) */
    lpcc->auto_corr[0] *= (1.0 + regular_term);

    /* Perform recursive calculation */
    if (LPC_LevinsonDurbinRecursion(lpcc, lpcc->auto_corr, coef_order, lpcc->parcor_coef, lpcc->error_vars) != LPC_ERROR_OK) {
        return LPC_ERROR_NG;
    }

    /* Consider the effect of windows in error variance */
    {
        uint32_t i;
        const double inverse_sqr = LPC_ComputeWindowInverseSquaredSum(window_type, num_samples);
        for (i = 0; i < coef_order + 1; i++) {
            lpcc->error_vars[i] *= inverse_sqr;
        }
    }

    return LPC_ERROR_OK;
}

/* Calculate LPC coefficients using Levinson-Durbin recursion (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficients(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *lpc_coef, uint32_t coef_order,
    LPCWindowType window_type, double regular_term)
{
    /* Argument check */
    if ((data == NULL) || (lpc_coef == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Check number of input samples */
    if (num_samples > lpcc->max_num_buffer_samples) {
        return LPC_APIRESULT_EXCEED_MAX_NUM_SAMPLES;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoef(lpcc, data, num_samples, coef_order, window_type, regular_term) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    /* If the calculation is successful, copy the result */
    memmove(lpc_coef, &lpcc->a_vecs[coef_order - 1][1], sizeof(double) * coef_order);

    return LPC_APIRESULT_OK;
}

/* Calculate all LPC coefficients up to a given order using Levinson-Durbin recursion (double precision) */
LPCApiResult LPCCalculator_CalculateMultipleLPCCoefficients(
    struct LPCCalculator* lpcc,
    const double* data, uint32_t num_samples, double **lpc_coefs, double *error_vars, uint32_t max_coef_order,
    LPCWindowType window_type, double regular_term)
{
    uint32_t k;

    /* Argument check */
    if ((data == NULL) || (lpc_coefs == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (max_coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Check number of input samples */
    if (num_samples > lpcc->max_num_buffer_samples) {
        return LPC_APIRESULT_EXCEED_MAX_NUM_SAMPLES;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoef(lpcc, data, num_samples, max_coef_order, window_type, regular_term) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    /* If the calculation is successful, copy the result */
    for (k = 0; k < max_coef_order; k++) {
        memmove(lpc_coefs[k], &lpcc->a_vecs[k][1], sizeof(double) * max_coef_order);
    }
    /* If the calculation is successful, copy the result */
    memmove(error_vars, lpcc->error_vars, sizeof(double) * (max_coef_order + 1));

    return LPC_APIRESULT_OK;
}

/* Cholesky decomposition */
static LPCError LPC_CholeskyDecomposition(
    double **Amat, int32_t dim, double *inv_diag)
{
    int32_t i, j, k;
    double sum;

    /* Argument check */
    assert((Amat != NULL) && (inv_diag != NULL));

    for (i = 0; i < dim; i++) {
        sum = Amat[i][i];
        for (k = i - 1; k >= 0; k--) {
            sum -= Amat[i][k] * Amat[i][k];
        }
        if (sum <= 0.0) {
            return LPC_ERROR_SINGULAR_MATRIX;
        }
        /* 1.0 / sqrt(sum) uses pow because division causes loss of digits */
        inv_diag[i] = pow(sum, -0.5);
        for (j = i + 1; j < dim; j++) {
            sum = Amat[i][j];
            for (k = i - 1; k >= 0; k--) {
                sum -= Amat[i][k] * Amat[j][k];
            }
            Amat[j][i] = sum * inv_diag[i];
        }
    }

    return LPC_ERROR_OK;
}

/* Solve Amat * xvec = bvec by Cholesky decomposition */
static LPCError LPC_SolveByCholeskyDecomposition(
        const double * const* Amat, int32_t dim, double *xvec, const double *bvec, const double *inv_diag)
{
    int32_t i, j;
    double sum;

    /* Argument check */
    assert((Amat != NULL) && (inv_diag != NULL) && (bvec != NULL) && (xvec != NULL));

    /* Solve linear equations using decomposition */
    for (i = 0; i < dim; i++) {
        sum = bvec[i];
        for (j = i - 1; j >= 0; j--) {
            sum -= Amat[i][j] * xvec[j];
        }
        xvec[i] = sum * inv_diag[i];
    }
    for (i = dim - 1; i >= 0; i--) {
        sum = xvec[i];
        for (j = i + 1; j < dim; j++) {
            sum -= Amat[j][i] * xvec[j];
        }
        xvec[i] = sum * inv_diag[i];
    }

    return LPC_ERROR_OK;
}

#if 1
/* Calculate the coefficient matrix using the auxiliary function method (forward residual) */
static LPCError LPCAF_CalculateCoefMatrixAndVector(
        const double *data, uint32_t num_samples,
        const double *a_vec, double **r_mat, double *r_vec,
        uint32_t coef_order, double *pobj_value)
{
    double obj_value;
    uint32_t smpl, i, j;

    assert(data != NULL);
    assert(a_vec != NULL);
    assert(r_mat != NULL);
    assert(r_vec != NULL);
    assert(pobj_value != NULL);
    assert(num_samples > coef_order);

    /* Initialize the matrix to 0 */
    for (i = 0; i < coef_order; i++) {
        r_vec[i] = 0.0;
        for (j = 0; j < coef_order; j++) {
            r_mat[i][j] = 0.0;
        }
    }

    obj_value = 0.0;

    for (smpl = coef_order; smpl < num_samples; smpl++) {
        /* Residual calculation */
        double residual = data[smpl];
        double inv_residual;
        for (i = 0; i < coef_order; i++) {
            residual += a_vec[i] * data[smpl - i - 1];
        }
        residual = fabs(residual);
        obj_value += residual;
        /* Round off residuals that are too small (avoid zero-ERO division, regularization) */
        residual = (residual < LPCAF_RESIDUAL_EPSILON) ? LPCAF_RESIDUAL_EPSILON : residual;
        inv_residual = 1.0 / residual;
        /* Accumulate into coefficient matrix */
        for (i = 0; i < coef_order; i++) {
            r_vec[i] -= data[smpl] * data[smpl - i - 1] * inv_residual;
            for (j = i; j < coef_order; j++) {
                r_mat[i][j] += data[smpl - i - 1] * data[smpl - j - 1] * inv_residual;
            }
        }
    }

    /* Extend to symmetric elements */
    for (i = 0; i < coef_order; i++) {
        for (j = i + 1; j < coef_order; j++) {
            r_mat[j][i] = r_mat[i][j];
        }
    }

    /* Set of objective function values ​​*/
    (*pobj_value) = obj_value / (num_samples - coef_order);

    return LPC_ERROR_OK;
}
#else
/* Calculate coefficient matrix using auxiliary function method (forward and backward residual) */
static LPCError LPCAF_CalculateCoefMatrixAndVector(
        const double *data, uint32_t num_samples,
        const double *a_vec, double **r_mat, double *r_vec,
        uint32_t coef_order, double *pobj_value)
{
    double obj_value;
    uint32_t smpl, i, j;

    assert(data != NULL);
    assert(a_vec != NULL);
    assert(r_mat != NULL);
    assert(r_vec != NULL);
    assert(pobj_value != NULL);
    assert(num_samples > coef_order);

    /* Initialize the matrix to 0 */
    for (i = 0; i < coef_order; i++) {
        r_vec[i] = 0.0;
        for (j = 0; j < coef_order; j++) {
            r_mat[i][j] = 0.0;
        }
    }

    obj_value = 0.0;

    for (smpl = coef_order; smpl < num_samples - coef_order; smpl++) {
        /* Residual calculation */
        double forward = data[smpl], backward = data[smpl];
        double inv_forward, inv_backward;
        for (i = 0; i < coef_order; i++) {
            forward += a_vec[i] * data[smpl - i - 1];
            backward += a_vec[i] * data[smpl + i + 1];
        }
        forward = fabs(forward);
        backward = fabs(backward);
        obj_value += (forward + backward);
        /* Round off residuals that are too small (avoid zero-ERO division, regularization) */
        forward = (forward < LPCAF_RESIDUAL_EPSILON) ? LPCAF_RESIDUAL_EPSILON : forward;
        backward = (backward < LPCAF_RESIDUAL_EPSILON) ? LPCAF_RESIDUAL_EPSILON : backward;
        inv_forward = 1.0 / forward;
        inv_backward = 1.0 / backward;
        /* Add to coefficient matrix */
        for (i = 0; i < coef_order; i++) {
            r_vec[i] -= data[smpl] * data[smpl - i - 1] * inv_forward;
            r_vec[i] -= data[smpl] * data[smpl + i + 1] * inv_backward;
            for (j = i; j < coef_order; j++) {
                r_mat[i][j] += data[smpl - i - 1] * data[smpl - j - 1] * inv_forward;
                r_mat[i][j] += data[smpl + i + 1] * data[smpl + j + 1] * inv_backward;
            }
        }
    }

    /* Extend to symmetric elements */
    for (i = 0; i < coef_order; i++) {
        for (j = i + 1; j < coef_order; j++) {
            r_mat[j][i] = r_mat[i][j];
        }
    }

    (*pobj_value) = obj_value / (2.0 * (num_samples - (2.0 * coef_order)));

    return LPC_ERROR_OK;
}
#endif

/* Coefficient calculation using auxiliary function method */
static LPCError LPC_CalculateCoefAF(
        struct LPCCalculator *lpcc, const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
        const uint32_t max_num_iteration, const double obj_epsilon, LPCWindowType window_type, double regular_term)
{
    uint32_t itr, i;
    double *a_vec = lpcc->work_buffer;
    double *r_vec = lpcc->u_vec;
    double **r_mat = lpcc->r_mat;
    double obj_value, prev_obj_value;
    LPCError err;

    /* Initialize coefficients using the Lebinson-Durbin method */
    if ((err = LPC_CalculateCoef(lpcc, data, num_samples, coef_order, window_type, regular_term)) != LPC_ERROR_OK) {
        return err;
    }
    memcpy(coef, &lpcc->a_vecs[coef_order - 1][1], sizeof(double) * coef_order);

    /*
/* If the zeroth-order autocorrelation (sum of squares of the signal) is small
* => Predict a silent output system with all coefficients set to 0 */
*/
    if (fabs(lpcc->auto_corr[0]) < FLT_EPSILON) {
        for (i = 0; i < coef_order + 1; i++) {
            lpcc->a_vecs[coef_order - 1][i] = 0.0;
        }
        return LPC_ERROR_OK;
    }

    prev_obj_value = FLT_MAX;
    for (itr = 0; itr < max_num_iteration; itr++) {
        /* 係数行列要素の計算 */
        if ((err = LPCAF_CalculateCoefMatrixAndVector(
                data, num_samples, a_vec, r_mat, r_vec, coef_order, &obj_value)) != LPC_ERROR_OK) {
            return err;
        }
        /* Cholesky decomposition */
        if ((err = LPC_CholeskyDecomposition(
                r_mat, (int32_t)coef_order, lpcc->v_vec)) == LPC_ERROR_SINGULAR_MATRIX) {
            /* In theory, a matrix becomes singular when all inputs are 0. Clear the coefficients to 0 and then exit. */
            for (i = 0; i < coef_order; i++) {
                lpcc->a_vecs[coef_order - 1][i] = 0.0;
            }
            return LPC_ERROR_OK;
        }
        /* Solve r_mat @ avec = r_vec using Cholesky decomposition */
        if ((err = LPC_SolveByCholeskyDecomposition(
                (const double * const *)r_mat, (int32_t)coef_order, coef, r_vec, lpcc->v_vec)) != LPC_ERROR_OK) {
            return err;
        }
        assert(err == LPC_ERROR_OK);
        /* Check convergence */
        if (fabs(prev_obj_value - obj_value) < obj_epsilon) {
            break;
        }
        prev_obj_value = obj_value;
    }

    return LPC_ERROR_OK;
}

/* Calculate LPC coefficients using auxiliary function method (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsAF(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    uint32_t max_num_iteration, LPCWindowType window_type, double regular_term)
{
    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (coef == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoefAF(lpcc, data, num_samples, coef, coef_order,
            max_num_iteration, 1e-8, window_type, regular_term) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    return LPC_APIRESULT_OK;
}

/* Coefficient calculation using Burg method */
static LPCError LPC_CalculateCoefBurg(
        struct LPCCalculator *lpcc, const double *data, uint32_t num_samples, double *coef, uint32_t coef_order)
{
#if 1
    uint32_t i, j, k;
    double *a_vec = lpcc->a_vecs[0];
    double **cov = lpcc->r_mat;
    LPCError err;

    /* Autocovariance matrix calculation */
    for (i = 0; i <= coef_order; i++) {
        if ((err = LPC_CalculateAutoCorrelation(
                        data, num_samples - i, &cov[i][i], coef_order + 1 - i)) != LPC_ERROR_OK) {
            return err;
        }
        for (j = i + 1; j <= coef_order; j++) {
            cov[j][i] = cov[i][j];
        }
    }

    /* Coefficient initialization */
    for (i = 0; i <= coef_order; i++) {
        a_vec[i] = 0.0;
    }
    a_vec[0] = 1.0;

    /* Calculate for each degree */
    for (k = 0; k < coef_order; k++) {
        double mu;
        double FkpBk = 0.0, sum = 0.0, Ck = 0.0;
        /* Fk + Bk */
        for (i = 0; i <= k; i++) {
            FkpBk += a_vec[i] * a_vec[i] * (cov[i][i] + cov[k + 1 - i][k + 1 - i]);
            /* Calculate only half of the non-diagonal elements using symmetry */
            for (j = i + 1; j <= k; j++) {
                sum += a_vec[i] * a_vec[j] * (cov[i][j] + cov[k + 1 - i][k + 1 - j]);
            }
        }
        FkpBk += 2.0 * sum;
        /* Ck */
        for (i = 0; i <= k; i++) {
            for (j = 0; j <= k; j++) {
                Ck += a_vec[i] * a_vec[j] * cov[i][k + 1 - j];
            }
        }
        /* Negative sign of reflection coefficient */
        mu = - 2.0 * Ck / FkpBk;
        assert(fabs(mu) <= 1.0);
        /* Update coefficients */
        for (i = 0; i <= (k + 1) / 2; i++) {
            double tmp1, tmp2;
            tmp1 = a_vec[i]; tmp2 = a_vec[k + 1 - i];
            a_vec[i]         = tmp1 + mu * tmp2;
            a_vec[k + 1 - i] = mu * tmp1 + tmp2;
        }
    }

    /* Set the solution */
    memcpy(coef, &lpcc->a_vecs[coef_order - 1][1], sizeof(double) * coef_order);
#else
    uint32_t i, k;
    double *a_vec = lpcc->a_vecs[0];
    double *f_vec, *b_vec;
    double Dk, mu;
    double tmp1, tmp2;

    /* Vector space allocation */
    f_vec = malloc(sizeof(double) * num_samples);
    b_vec = malloc(sizeof(double) * num_samples);

    /* Initialize each vector */
    for (k = 0; k < coef_order + 1; k++) {
        a_vec[k] = 0.0;
    }
    a_vec[0] = 1.0;
    memcpy(f_vec, data, sizeof(double) * num_samples);
    memcpy(b_vec, data, sizeof(double) * num_samples);

    /* Initialize Dk */
    Dk = 0.0;
    for (i = 0; i < num_samples; i++) {
        Dk += 2.0 * f_vec[i] * f_vec[i];
    }
    Dk -= f_vec[0] * f_vec[0] + f_vec[num_samples - 1] * f_vec[num_samples - 1];

    /* Burg recursion algorithm */
    for (k = 0; k < coef_order; k++) {
        /* Calculate the reflection (PARCOR) coefficient */
        mu = 0.0;
        for (i = 0; i < num_samples - k - 1; i++) {
            mu += f_vec[i + k + 1] * b_vec[i];
        }
        mu *= -2.0 / Dk;
        assert(fabs(mu) < 1.0);
        /* Update a_vec */
        for (i = 0; i <= (k + 1) / 2; i++) {
            tmp1 = a_vec[i]; tmp2 = a_vec[k + 1 - i];
            a_vec[i]         = tmp1 + mu * tmp2;
            a_vec[k + 1 - i] = mu * tmp1 + tmp2;
        }
        /* Update f_vec, b_vec */
        for (i = 0; i < num_samples - k - 1; i++) {
            tmp1 = f_vec[i + k + 1]; tmp2 = b_vec[i];
            f_vec[i + k + 1] = tmp1 + mu * tmp2;
            b_vec[i]         = mu * tmp1 + tmp2;
        }
        /* Update Dk */
        Dk = (1.0 - mu * mu) * Dk - f_vec[k + 1] * f_vec[k + 1] - b_vec[num_samples - k - 2] * b_vec[num_samples - k - 2];
    }

    /* Copy coefficients */
    memcpy(a_vec, &lpcc->a_vecs[coef_order - 1][1], sizeof(double) * coef_order);

    free(b_vec);
    free(f_vec);
#endif
    return LPC_ERROR_OK;
}

/* Calculate LPC coefficients using Burg method (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsBurg(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order)
{
    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (coef == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoefBurg(lpcc, data, num_samples, coef, coef_order) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    return LPC_APIRESULT_OK;
}

/* Calculate the covariance matrix */
static LPCError LPCSVR_CalculateCovarianceMatrix(
    const double *data, uint32_t num_samples, double **cov, uint32_t dim)
{
    uint32_t i, j, smpl;

    /* Argument check */
    if ((data == NULL) || (cov == NULL)) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    for (i = 0; i < dim; i++) {
        for (j = i; j < dim; j++) {
            cov[i][j] = 0.0;
        }
    }

    for (smpl = 0; smpl < num_samples - dim; smpl++) {
        const double *pdata = &data[smpl];
        for (i = 0; i < dim; i++) {
            const double s = pdata[i];
            for (j = i; j < dim; j++) {
                cov[i][j] += s * pdata[j];
            }
        }
    }

    for (i = 0; i < dim; i++) {
        for (j = i + 1; j < dim; j++) {
            cov[j][i] = cov[i][j];
        }
    }

    return LPC_ERROR_OK;
}

/* Average code length of Recursive Golomb-Rice code */
static double LPCSVR_CalculateRGRMeanCodeLength(double mean_abs_error, uint32_t bps)
{
    const double intmean = mean_abs_error * (1 << bps); /* Average value when quantized to an integer */
    const double rho = 1.0 / (1.0 + intmean);
    const uint32_t k2 = (uint32_t)LPC_MAX(0, LPC_Log2(log(0.5127629514) / log(1.0 - rho)));
    const uint32_t k1 = k2 + 1;
    const double k1factor = pow(1.0 - rho, (double)(1 << k1));
    const double k2factor = pow(1.0 - rho, (double)(1 << k2));
    return (1.0 + k1) * (1.0 - k1factor) + (1.0 + k2 + (1.0 / (1.0 - k2factor))) * k1factor;
}

/* Coefficient calculation by SVR */
static LPCError LPC_CalculateCoefSVR(
    struct LPCCalculator *lpcc, const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    const uint32_t max_num_iteration, const double obj_epsilon, LPCWindowType window_type,
    double regular_term, const double *margin_list, uint32_t margin_list_size)
{
#define BITS_PER_SAMPLE 16
    uint32_t itr, i, j, smpl;
    double *r_vec = lpcc->u_vec;
    double *low = lpcc->v_vec;
    double *best_coef = lpcc->work_buffer;
    double *delta = lpcc->parcor_coef;
    double *init_coef = lpcc->auto_corr;
    double **cov = lpcc->r_mat;
    double *residual = lpcc->buffer;
    double obj_value, prev_obj_value, min_obj_value;
    LPCError err;

    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (margin_list == NULL) || (margin_list_size == 0)) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /* If not learning, exit immediately */
    if (max_num_iteration == 0) {
        return LPC_ERROR_OK;
    }

    /* Calculate the covariance matrix */
    if ((err = LPCSVR_CalculateCovarianceMatrix(data, num_samples, cov, coef_order)) != LPC_ERROR_OK) {
        return err;
    }
    /* Ridge regularization */
    for (i = 0; i < coef_order; i++) {
        cov[i][i] *= (1.0 + regular_term);
    }
    /* Cholesky decomposition */
    if ((err = LPC_CholeskyDecomposition(cov, (int32_t)coef_order, low)) == LPC_ERROR_SINGULAR_MATRIX) {
        /* In theory, a matrix becomes singular when all inputs are 0. Clear the coefficients to 0 and then exit. */
        for (i = 0; i < coef_order; i++) {
            coef[i] = 0.0;
        }
        return LPC_ERROR_OK;
    }

    /* Set the initial value to the argument coefficient */
    memcpy(init_coef, coef, sizeof(double) * coef_order);
    memcpy(best_coef, init_coef, sizeof(double) * coef_order);

    /* TODO: It seems that the calculation of residuals will be faster if the order of coefficients is reversed (needs verification) */

    min_obj_value = FLT_MAX;
    for (j = 0; j < margin_list_size; j++) {
        const double margin = margin_list[j];
        prev_obj_value = FLT_MAX;
        memcpy(coef, init_coef, sizeof(double) * coef_order);
        for (itr = 0; itr < max_num_iteration; itr++) {
            double mabse = 0.0;
            /* Residual calculation/residual soft threshold */
            memcpy(residual, data, sizeof(double) * num_samples);
            for (i = 0; i < coef_order; i++) {
                r_vec[i] = 0.0;
            }
            for (smpl = coef_order; smpl < num_samples; smpl++) {
                for (i = 0; i < coef_order; i++) {
                    residual[smpl] += coef[i] * data[smpl - i - 1];
                }
                mabse += LPC_ABS(residual[smpl]);
                residual[smpl] = LPC_SOFT_THRESHOLD(residual[smpl], margin);
                for (i = 0; i < coef_order; i++) {
                    r_vec[i] += residual[smpl] * data[smpl - i - 1];
                }
            }
            obj_value = LPCSVR_CalculateRGRMeanCodeLength(mabse / num_samples, BITS_PER_SAMPLE);
            /* Solve cov @ delta = r_vec using Cholesky decomposition */
            if ((err = LPC_SolveByCholeskyDecomposition(
                    (const double * const *)cov, (int32_t)coef_order, delta, r_vec, low)) != LPC_ERROR_OK) {
                return err;
            }
            /* Update best coefficients */
            if (obj_value < min_obj_value) {
                memcpy(best_coef, coef, sizeof(double) * coef_order);
                min_obj_value = obj_value;
            }
            /* Check convergence */
            if ((prev_obj_value < obj_value) || (fabs(prev_obj_value - obj_value) < obj_epsilon)) {
                break;
            }
            /* Update coefficients */
            for (i = 0; i < coef_order; i++) {
                coef[i] += delta[i];
            }
            prev_obj_value = obj_value;
        }
    }

    /* Record the best coefficient */
    memcpy(coef, best_coef, sizeof(double) * coef_order);

    return LPC_ERROR_OK;
#undef BITS_PER_SAMPLE
}

/* Calculate LPC coefficients from SVR (double precision) */
LPCApiResult LPCCalculator_CalculateLPCCoefficientsSVR(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, double *coef, uint32_t coef_order,
    uint32_t max_num_iteration, LPCWindowType window_type,
    double regular_term, const double *margin_list, uint32_t margin_list_size)
{
    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (coef == NULL)
        || (margin_list == NULL) || (margin_list_size == 0)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoefSVR(lpcc, data, num_samples, coef, coef_order,
            max_num_iteration, 1e-8, window_type, regular_term, margin_list, margin_list_size) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    return LPC_APIRESULT_OK;
}

/* Calculate the estimated code length per sample from the input data */
LPCApiResult LPCCalculator_EstimateCodeLength(
        struct LPCCalculator *lpcc,
        const double *data, uint32_t num_samples, uint32_t bits_per_sample,
        uint32_t coef_order, double *length_per_sample_bits, LPCWindowType window_type)
{
    uint32_t ord;
    double log2_mean_res_power, log2_var_ratio;

    /* constant value */
#define BETA_CONST_FOR_LAPLACE_DIST   (1.9426950408889634)  /* sqrt(2 * E * E) */
#define BETA_CONST_FOR_GAUSS_DIST     (2.047095585180641)   /* sqrt(2 * E * PI) */

    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (length_per_sample_bits == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoef(lpcc, data, num_samples, coef_order, window_type, 0.0) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    /* Calculate log2(power average) */
    log2_mean_res_power = lpcc->auto_corr[0]; /* 0th order sample autocorrelation is power */
    /* Convert to integer PCM amplitude (guarantee double density) */
    log2_mean_res_power *= pow(2, (double)(2.0 * (bits_per_sample - 1)));
    if (fabs(log2_mean_res_power) <= FLT_MIN) {
        /* If there is almost no sound, set the code length to 0 */
        (*length_per_sample_bits) = 0.0;
        return LPC_APIRESULT_OK;
    }
    log2_mean_res_power = LPC_Log2((double)log2_mean_res_power) - LPC_Log2((double)num_samples);

    /* Calculate sum(log2(1 - (parcor * parcor))) */
    log2_var_ratio = 0.0;
    for (ord = 0; ord < coef_order; ord++) {
        log2_var_ratio += LPC_Log2(1.0 - lpcc->parcor_coef[ord] * lpcc->parcor_coef[ord]);
    }

    /* Entropy calculation */
    /* → This gives the minimum number of bits per sample */
    (*length_per_sample_bits) = BETA_CONST_FOR_LAPLACE_DIST + 0.5f * (log2_mean_res_power + log2_var_ratio);

    /* If the estimated number of bits is negative, we hope to be able to encode with 1 bit per sample */
    /* Note) In this case the input audio power is very low */
    if ((*length_per_sample_bits) <= 0) {
        (*length_per_sample_bits) = 1.0;
        return LPC_APIRESULT_OK;
    }

#undef BETA_CONST_FOR_LAPLACE_DIST
#undef BETA_CONST_FOR_GAUSS_DIST

    return LPC_APIRESULT_OK;
}

/* Calculate the MDL (Minimum Description Length) */
LPCApiResult LPCCalculator_CalculateMDL(
    struct LPCCalculator *lpcc,
    const double *data, uint32_t num_samples, uint32_t coef_order, double *mdl,
    LPCWindowType window_type)
{
    uint32_t k;
    double tmp;

    /* Argument check */
    if ((lpcc == NULL) || (data == NULL) || (mdl == NULL)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Coefficient calculation */
    if (LPC_CalculateCoef(lpcc, data, num_samples, coef_order, window_type, 0.0) != LPC_ERROR_OK) {
        return LPC_APIRESULT_FAILED_TO_CALCULATION;
    }

    /* Calculate the first term */
    /* The first coefficient is definitely 0, so skip it */
    tmp = 0.0;
    for (k = 1; k <= coef_order; k++) {
        tmp += log(1.0 - lpcc->parcor_coef[k] * lpcc->parcor_coef[k]);
    }
    tmp *= num_samples;

    /* Calculate the second term */
    tmp += coef_order * log(num_samples);

    (*mdl) = tmp;

    return LPC_APIRESULT_OK;
}

/* Convert LPC coefficients to PARCOR coefficients */
static LPCError LPC_ConvertLPCtoPARCORDouble(
    struct LPCCalculator *lpcc, const double *lpc_coef, uint32_t coef_order, double *parcor_coef)
{
    int32_t i, k;
    double *tmplpc_coef, *a_vec;

    /* Argument check */
    if ((lpcc == NULL) || (lpc_coef == NULL) || (parcor_coef == NULL)) {
        return LPC_ERROR_INVALID_ARGUMENT;
    }

    /* Check degree */
    assert(coef_order <= lpcc->max_order);

    /* Allocate work area */
    tmplpc_coef = lpcc->work_buffer;
    a_vec = lpcc->a_vecs[0];

    memcpy(tmplpc_coef, lpc_coef, sizeof(double) * coef_order);

    /* Convert to PARCOR coefficients */
    for (i = (int32_t)(coef_order - 1); i >= 0; i--) {
        const double gamma = tmplpc_coef[i];
        assert(fabs(gamma) < 1.0);
        parcor_coef[i] = -gamma;
        for (k = 0; k < i; k++) {
            a_vec[k] = tmplpc_coef[k];
        }
        for (k = 0; k < i; k++) {
            tmplpc_coef[k] = (a_vec[k] - gamma * a_vec[i - k - 1]) / (1.0 - gamma * gamma);
        }
    }

    return LPC_ERROR_OK;
}

/* Convert LPC coefficients to PARCOR coefficients and quantize */
LPCApiResult LPC_QuantizeCoefficientsAsPARCOR(
    struct LPCCalculator *lpcc,
    const double *lpc_coef, uint32_t coef_order, uint32_t nbits_precision, int32_t *int_coef)
{
    uint32_t ord;
    int32_t qtmp;
    const int32_t qmax = (1 << (nbits_precision - 1));

    /* Argument check */
    if ((lpcc == NULL) || (lpc_coef == NULL)
        || (int_coef == NULL) || (nbits_precision == 0)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Check degree */
    if (coef_order > lpcc->max_order) {
        return LPC_APIRESULT_EXCEED_MAX_ORDER;
    }

    /* Convert to PARCOR coefficients */
    if (LPC_ConvertLPCtoPARCORDouble(lpcc, lpc_coef, coef_order, lpcc->parcor_coef) != LPC_ERROR_OK) {
        return LPC_APIRESULT_NG;
    }

    /* Quantize and output PARCOR coefficients */
    for (ord = 0; ord < coef_order; ord++) {
        assert(fabs(lpcc->parcor_coef[ord]) < 1.0);
        qtmp = (int32_t)LPC_Round(lpcc->parcor_coef[ord] * pow(2.0, nbits_precision - 1));
        /* Rounding to positive and negative boundaries */
        if (qtmp >= qmax) {
            qtmp = qmax - 1;
        } else if (qtmp < -qmax) {
            qtmp = -qmax;
        }
        int_coef[ord] = qtmp;
    }

    return LPC_APIRESULT_OK;
}

/* Integer quantization of LPC coefficients */
LPCApiResult LPC_QuantizeCoefficients(
    const double *double_coef, uint32_t coef_order, uint32_t nbits_precision, uint32_t max_bits,
    int32_t *int_coef, uint32_t *coef_rshift)
{
    uint32_t rshift;
    int32_t ord, ndigit, qtmp;
    double max, qerror;
    const int32_t qmax = (1 << (nbits_precision - 1));

    /* Argument check */
    if ((double_coef == NULL) || (int_coef == NULL)
            || (coef_rshift == NULL) || (nbits_precision == 0)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Calculate the absolute value of the coefficient */
    max = 0.0;
    for (ord = 0; ord < (int32_t)coef_order; ord++) {
        if (max < fabs(double_coef[ord])) {
            max = fabs(double_coef[ord]);
        }
    }

    /* If the value is too small to be represented in the given number of bits, it is considered to be 0 */
    if (max <= pow(2.0, -(int32_t)(nbits_precision - 1))) {
        (*coef_rshift) = nbits_precision;
        memset(int_coef, 0, sizeof(int32_t) * coef_order);
        return LPC_APIRESULT_OK;
    }

    /* Calculate the right shift amount to put the maximum value in [1/2, 1) */
    /* Calculate max = x * 2^ndigit, |x| in [1/2, 1) */
    (void)frexp(max, &ndigit);
    /* Drop the sign bit */
    nbits_precision--;
    /* Calculate the shift amount to make it representable with nbits_precision */
    assert((int32_t)nbits_precision >= ndigit);
    rshift = (uint32_t)((int32_t)nbits_precision - ndigit);

    /* If the right shift amount exceeds the maximum, truncate (if the maximum coefficient value is small) */
    if (rshift >= max_bits) {
        rshift = max_bits - 1;
    }

    /* Quantization */
    qerror = 0.0;
    for (ord = (int32_t)coef_order - 1; ord >= 0; ord--) {
        /* Incorporate the error of the previous coefficient and quantize it */
        /* We don't want to introduce errors at the beginning of the impulse, so we process from the end */
        qerror += double_coef[ord] * pow(2.0, rshift);
        qtmp = (int32_t)LPC_Round(qerror);
        /* Rounding to positive and negative boundaries */
        if (qtmp >= qmax) {
            qtmp = qmax - 1;
        } else if (qtmp < -qmax) {
            qtmp = -qmax;
        }
        /* The subtracted amount remains as the quantization error */
        qerror -= qtmp;
        int_coef[ord] = qtmp;
    }
    (*coef_rshift) = rshift;

    return LPC_APIRESULT_OK;
}

/* Prediction/error output by LPC coefficients */
LPCApiResult LPC_Predict(
    const int32_t *data, uint32_t num_samples,
    const int32_t *coef, uint32_t coef_order, int32_t *residual, uint32_t coef_rshift)
{
    uint32_t smpl, ord;

    /* Argument check */
    if ((data == NULL) || (coef == NULL)
            || (residual == NULL) || (coef_rshift == 0)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    memcpy(residual, data, sizeof(int32_t) * num_samples);

    /* Prediction based on LPC coefficients */
    for (smpl = 1; smpl < coef_order; smpl++) {
        int32_t predict = (1 << (coef_rshift - 1));
        for (ord = 0; ord < smpl; ord++) {
            predict += (coef[ord] * data[smpl - ord - 1]);
        }
        residual[smpl] += (predict >> coef_rshift);
    }
    for (smpl = coef_order; smpl < num_samples; smpl++) {
        int32_t predict = (1 << (coef_rshift - 1));
        for (ord = 0; ord < coef_order; ord++) {
            predict += (coef[ord] * data[smpl - ord - 1]);
        }
        residual[smpl] += (predict >> coef_rshift);
    }

    return LPC_APIRESULT_OK;
}

/* Synthesis by LPC coefficients (in-place) */
LPCApiResult LPC_Synthesize(
    int32_t *data, uint32_t num_samples,
    const int32_t *coef, uint32_t coef_order, uint32_t coef_rshift)
{
    uint32_t smpl, ord;

    /* Argument check */
    if ((data == NULL) || (coef == NULL) || (coef_rshift == 0)) {
        return LPC_APIRESULT_INVALID_ARGUMENT;
    }

    /* Prediction based on LPC coefficients */
    for (smpl = 1; smpl < coef_order; smpl++) {
        int32_t predict = (1 << (coef_rshift - 1));
        for (ord = 0; ord < smpl; ord++) {
            predict += (coef[ord] * data[smpl - ord - 1]);
        }
        data[smpl] -= (predict >> coef_rshift);
    }
    for (smpl = coef_order; smpl < num_samples; smpl++) {
        int32_t predict = (1 << (coef_rshift - 1));
        for (ord = 0; ord < coef_order; ord++) {
            predict += (coef[ord] * data[smpl - ord - 1]);
        }
        data[smpl] -= (predict >> coef_rshift);
    }

    return LPC_APIRESULT_OK;
}
