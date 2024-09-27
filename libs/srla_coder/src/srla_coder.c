#include "srla_coder.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "srla_internal.h"
#include "srla_utility.h"

/* Use macro expansion */
#define SRLACODER_USE_MACROS 1

/* memory alignment */
#define SRLACODER_MEMORY_ALIGNMENT 16
/* log2(maximum number of divisions) */
#define SRLACODER_LOG2_MAX_NUM_PARTITIONS 10
/* Maximum number of divisions */
#define SRLACODER_MAX_NUM_PARTITIONS (1 << SRLACODER_LOG2_MAX_NUM_PARTITIONS)
/* Number of bits in parameter recording area */
#define SRLACODER_RICE_PARAMETER_BITS 5
/* Gamma code length size */
#define SRLACODER_GAMMA_BITS(uint) (((uint) == 0) ? 1 : ((2 * SRLAUTILITY_LOG2CEIL(uint + 2)) - 1))

/* Sign Sensitivity */
typedef enum SRLACoderCodeTypeTag {
    SRLACODER_CODE_TYPE_RICE = 0, /* TODO: Maybe make it Golomb in the future */
    SRLACODER_CODE_TYPE_RECURSIVE_RICE = 1,
    SRLACODER_CODE_TYPE_ALLZERO = 2, /* In the case of pseudo-stereo/reverse phase, one channel may be all 0 */
    SRLACODER_CODE_TYPE_INVALID
} SRLACoderCodeType;

/* encoding handle */
struct SRLACoder {
    uint8_t alloced_by_own;
    double part_mean[SRLACODER_LOG2_MAX_NUM_PARTITIONS + 1][SRLACODER_MAX_NUM_PARTITIONS];
    uint32_t *uval_buffer;
    void *work;
};

/* Calculate the work size required to create the encoding handle */
int32_t SRLACoder_CalculateWorkSize(uint32_t max_num_samples)
{
    int32_t work_size;

    /* Size of handle */
    work_size = sizeof(struct SRLACoder) + SRLACODER_MEMORY_ALIGNMENT;

    /* Size of code sample buffer */
    work_size += SRLACODER_MEMORY_ALIGNMENT + sizeof(uint32_t) * max_num_samples;

    return work_size;
}

/* Create an encoding handle */
struct SRLACoder* SRLACoder_Create(uint32_t max_num_samples, void *work, int32_t work_size)
{
    struct SRLACoder *coder;
    uint8_t tmp_alloc_by_own = 0;
    uint8_t *work_ptr;

    /* In case of pre-allocation of work area */
    if ((work == NULL) && (work_size == 0)) {
        /* Replace arguments with your own calculated values ​​*/
        if ((work_size = SRLACoder_CalculateWorkSize(max_num_samples)) < 0) {
            return NULL;
        }
        work = malloc((uint32_t)work_size);
        tmp_alloc_by_own = 1;
    }

    /* Argument check */
    if ((work == NULL) || (work_size < SRLACoder_CalculateWorkSize(max_num_samples))) {
        return NULL;
    }

    /* Get the start of the work area */
    work_ptr = (uint8_t *)work;

    /* Allocate handle area */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLACODER_MEMORY_ALIGNMENT);
    coder = (struct SRLACoder *)work_ptr;
    work_ptr += sizeof(struct SRLACoder);

    /* Allocate code sample buffer */
    work_ptr = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLACODER_MEMORY_ALIGNMENT);
    coder->uval_buffer = (uint32_t *)work_ptr;
    work_ptr += sizeof(uint32_t) * max_num_samples;

    /* Handle member settings */
    coder->alloced_by_own = tmp_alloc_by_own;
    coder->work = work;

    /* Initialize division information */
    {
        uint32_t i, j;
        for (i = 0; i < SRLACODER_LOG2_MAX_NUM_PARTITIONS + 1; i++) {
            for (j = 0; j < SRLACODER_MAX_NUM_PARTITIONS; j++) {
                coder->part_mean[i][j] = UINT32_MAX;
            }
        }
    }

    return coder;
}

/* Destroy the encoding handle */
void SRLACoder_Destroy(struct SRLACoder *coder)
{
    if (coder != NULL) {
        /* Release the area if you have it yourself */
        if (coder->alloced_by_own == 1) {
            free(coder->work);
        }
    }
}

/* Output gamma code */
static void Gamma_PutCode(struct BitStream *stream, uint32_t val)
{
    uint32_t ndigit;

    SRLA_ASSERT(stream != NULL);

    if (val == 0) {
        /* If the encoding target is 0, output 1 and exit */
        BitWriter_PutBits(stream, 1, 1);
        return;
    }

    /* Get the number of digits */
    ndigit = SRLAUTILITY_LOG2CEIL(val + 2);
    /* Continue with 0 for number of digits -1 */
    BitWriter_PutBits(stream, 0, ndigit - 1);
    /* Output the codeword in binary using the number of digits */
    BitWriter_PutBits(stream, val + 1, ndigit);
}

/* Get gamma code */
static uint32_t Gamma_GetCode(struct BitStream *stream)
{
    uint32_t ndigit;
    uint32_t bitsbuf;

    SRLA_ASSERT(stream != NULL);

    /* Get the number of digits */
    /* Increase the number of digits until a 1 appears */
    BitReader_GetZeroRunLength(stream, &ndigit);
    /* At least 1, so we'll make it work */
    ndigit++;

    /* 0 if the number of digits is 1 */
    if (ndigit == 1) {
        return 0;
    }

    /* Output code word from digits */
    BitReader_GetBits(stream, &bitsbuf, ndigit - 1);
    return (uint32_t)((1UL << (ndigit - 1)) + bitsbuf - 1);
}

/* Output of Rice code */
static void Rice_PutCode(struct BitStream *stream, uint32_t k, uint32_t uval)
{
    SRLA_ASSERT(stream != NULL);

    BitWriter_PutZeroRun(stream, uval >> k);
    BitWriter_PutBits(stream, uval, k);
}

/* Output of recursive Rice code */
#if defined(SRLACODER_USE_MACROS)
#define RecursiveRice_PutCode(stream, k1, k2, uval)\
    do {\
        const uint32_t k1pow__ = 1U << (k1);\
        \
        SRLA_ASSERT((stream) != NULL);\
        \
        if ((uval) < k1pow__) {\
            /* Encoding in the first stage */\
            BitWriter_PutBits((stream), k1pow__ | (uval), (k1) + 1);\
        } else {\
            /* Subtract with the first stage parameters, and Rice encode with the second stage parameters */\
            const uint32_t tmp_uval__ = (uval) - k1pow__;\
            BitWriter_PutZeroRun((stream), 1 + (tmp_uval__ >> (k2)));\
            BitWriter_PutBits((stream), tmp_uval__, (k2));\
        }\
    } while (0);
#else
static void RecursiveRice_PutCode(struct BitStream *stream, uint32_t k1, uint32_t k2, uint32_t uval)
{
    const uint32_t k1pow = 1U << k1;

    SRLA_ASSERT(stream != NULL);

    if (uval < k1pow) {
        /* Encoding in the first stage */
        BitWriter_PutBits(stream, k1pow | uval, k1 + 1);
    } else {
        /* Subtract with the first stage parameters, and Rice encode with the second stage parameters */
        uval -= k1pow;
        BitWriter_PutZeroRun(stream, 1 + (uval >> k2));
        BitWriter_PutBits(stream, uval, k2);
    }
}
#endif

/* Get Rice code */
static uint32_t Rice_GetCode(struct BitStream *stream, uint32_t k)
{
    uint32_t quot, uval;

    SRLA_ASSERT(stream != NULL);

    /* Get the quotient (alpha sign) */
    BitReader_GetZeroRunLength(stream, &quot);

    /* Get the remainder */
    BitReader_GetBits(stream, &uval, k);

    return (quot << k) + uval;
}

/* Get recursive Rice code */
#if defined(SRLACODER_USE_MACROS)
#define RecursiveRice_GetCode(stream, k1, k2, uval)\
    do {\
        uint32_t quot__;\
        \
        SRLA_ASSERT((stream) != NULL);\
        SRLA_ASSERT((uval) != NULL);\
        SRLA_ASSERT((k1) == ((k2) + 1));\
        \
        /* Get the business section */\
        BitReader_GetZeroRunLength((stream), &quot__);\
        \
        /* Get the remainder */\
        BitReader_GetBits(stream, uval, (k2) + !(quot__));\
        (*uval) |= ((quot__ + !!(quot__)) << (k2));\
    } while (0);
#else
static void RecursiveRice_GetCode(struct BitStream *stream, uint32_t k1, uint32_t k2, uint32_t *uval)
{
    uint32_t quot;

    SRLA_ASSERT(stream != NULL);
    SRLA_ASSERT(uval != NULL);
    SRLA_ASSERT(k1 == (k2 + 1));

    /* Get the business section */
    BitReader_GetZeroRunLength(stream, &quot);

    /* Get the remainder */
    BitReader_GetBits(stream, uval, k2 + !(quot));
    (*uval) |= ((quot + !!(quot)) << (k2));
}
#endif

/* Calculate optimal encoding parameters */
static void SRLACoder_CalculateOptimalRiceParameter(
    const double mean, uint32_t *optk, double *bits_per_sample)
{
    uint32_t k;
    double rho, fk, bps;
#define OPTX 0.5127629514437670454896078808815218508243560791015625 /* Solution of (x - 1)^2 + ln(2) x ln(x) = 0 */

    /* Maximum likelihood estimation of geometric distribution parameters */
    rho = 1.0 / (1.0 + mean);

    /* Calculate the optimal parameters */
    k = (uint32_t)SRLAUTILITY_MAX(0, SRLAUtility_Round(SRLAUtility_Log2(log(OPTX) / log(1.0 - rho))));

    /* Calculate the average code length */
    fk = pow(1.0 - rho, (double)(1 << k));
    bps = k + (1.0 / (1.0 - fk));

    /* Result output */
    (*optk) = k;

    if (bits_per_sample != NULL) {
        (*bits_per_sample) = bps;
    }

#undef OPTX
}

/* Calculate the partial derivative with respect to k1 */
static double SRLACoder_CalculateMeanCodelength(double rho, uint32_t k1, uint32_t k2)
{
    const double fk1 = pow(1.0 - rho, (double)(1 << k1));
    const double fk2 = pow(1.0 - rho, (double)(1 << k2));
    return (1.0 + k1) * (1.0 - fk1) + (1.0 + k2 + (1.0 / (1.0 - fk2))) * fk1;
}

#if 0
/* Calculate the partial derivative with respect to k1 */
static double SRLACoder_CalculateDiffk1(double rho, double k1, double k2)
{
    const double k1pow = pow(2.0, k1);
    const double fk1 = pow(1.0 - rho, k1pow);
    const double fk2 = pow(1.0 - rho, pow(2.0, k2));
    const double fk1d = k1pow * fk1 * log(1.0 - rho) * log(2.0);
    return 1.0 - fk1 + (k2 - k1 + 1.0 / (1.0 - fk2)) * fk1d;
}
#endif

/* Calculate optimal encoding parameters */
static void SRLACoder_CalculateOptimalRecursiveRiceParameter(
    const double mean, uint32_t *optk1, uint32_t *optk2, double *bits_per_sample)
{
    uint32_t k1, k2;
    double rho;
#define OPTX 0.5127629514437670454896078808815218508243560791015625 /* Solution of (x - 1)^2 + ln(2) x ln(x) = 0 */

    /* Maximum likelihood estimation of geometric distribution parameters */
    rho = 1.0 / (1.0 + mean);

    /* Calculate the optimal parameters */
#if 0
    k2 = (uint32_t)SRLAUTILITY_MAX(0, floor(SRLAUtility_Log2(log(OPTX) / log(1.0 - rho))));
    k1 = k2 + 1;
#else
    /* Fast approximate calculation */
    {
#define MLNOPTX (0.66794162356) /* -ln(OPTX) */
        const uint32_t opt_golomb_param = (uint32_t)SRLAUTILITY_MAX(1, MLNOPTX * (1.0 + mean));
        k2 = SRLAUTILITY_LOG2FLOOR(opt_golomb_param);
        k1 = k2 + 1;
    }
#endif

#if 0
    /* Find k1 using bisection */
    /* Note: In general, k1 = k2 + 1 does not hold, which is disadvantageous for encoding! */
    {
        uint32_t i;
        double k1tmp, d1, d2;
        const double k2tmp = SRLAUtility_Log2(log(OPTX) / log(1.0 - rho));
        double k1min = k2tmp - 1.5;
        double k1max = k2tmp + 1.5;
        for (i = 0; i < 5; i++) {
            k1tmp = (k1max + k1min) / 2.0;
            d1 = SRLACoder_CalculateDiffk1(rho, k1tmp, k2tmp);
            d2 = SRLACoder_CalculateDiffk1(rho, k1min, k2tmp);
            if (SRLAUTILITY_SIGN(d1) == SRLAUTILITY_SIGN(d2)) {
                k1min = k1tmp;
            } else {
                k1max = k1tmp;
            }
        }
        k1 = (uint32_t)SRLAUTILITY_MAX(0, ceil(k1tmp));
        k2 = (uint32_t)SRLAUTILITY_MAX(0, floor(k2tmp));
    }
#endif

    /* Result output */
    (*optk1) = k1;
    (*optk2) = k2;

    /* Calculate the average code length */
    if (bits_per_sample != NULL) {
        (*bits_per_sample) = SRLACoder_CalculateMeanCodelength(rho, k1, k2);
    }

#undef OPTX
}

/* Output Rice code length */
static uint32_t Rice_GetCodeLength(uint32_t k, uint32_t uval)
{
    return 1 + k  + (uval >> k);
}

/* Calculate the recursive Rice code length for an array */
static uint32_t RecursiveRice_ComputeCodeLength(const uint32_t *data, uint32_t num_samples, uint32_t k1, uint32_t k2)
{
    uint32_t smpl, length;
    const uint32_t k1pow = 1U << k1;

    SRLA_ASSERT(data != NULL);

#if 0
    length = 0;
    for (smpl = 0; smpl < num_samples; smpl++) {
        const uint32_t uval = data[smpl];
        if (uval < k1pow) {
            /* Encoding in the first stage */
            length += (k1 + 1);
        } else {
            /* Subtract with the first stage parameters, and Rice encode with the second stage parameters */
            length += (k2 + 2 + ((uval - k1pow) >> k2));
        }
    }
#else
    SRLA_ASSERT((k2 + 1) == k1);
    length = (k1 + 1) * num_samples;
    for (smpl = 0; smpl < num_samples; smpl++) {
        length += (SRLAUTILITY_MAX(0, (int32_t)data[smpl] - (int32_t)k1pow) >> k2);
    }
#endif

    return length;
}

static void SRLACoder_SearchBestCodeAndPartition(
    struct SRLACoder *coder, const int32_t *data, uint32_t num_samples,
    SRLACoderCodeType *code_type, uint32_t *best_partition_order, uint32_t *best_code_length)
{
    uint32_t max_porder, max_num_partitions;
    uint32_t porder, part, best_porder, min_bits, smpl, max_uval;
    SRLACoderCodeType tmp_code_type = SRLACODER_CODE_TYPE_INVALID;

    /* Determine the maximum number of divisions */
    max_porder = 1;
    while ((num_samples % (1 << max_porder)) == 0) {
        max_porder++;
    }
    max_porder = SRLAUTILITY_MIN(max_porder - 1, SRLACODER_LOG2_MAX_NUM_PARTITIONS);
    max_num_partitions = (1 << max_porder);

    /* Calculate the average for each fold */
    {
        int32_t i;

        /* Average value at finest division */
        max_uval = 0;
        for (part = 0; part < max_num_partitions; part++) {
            const uint32_t nsmpl = num_samples / max_num_partitions;
            double part_sum = 0.0;
            for (smpl = 0; smpl < nsmpl; smpl++) {
                /* Cache the uint32 conversion result */
                const uint32_t uval = SRLAUTILITY_SINT32_TO_UINT32(data[part * nsmpl + smpl]);
                coder->uval_buffer[part * nsmpl + smpl] = uval;
                part_sum += uval;
                max_uval = SRLAUTILITY_MAX(max_uval, uval);
            }
            coder->part_mean[max_porder][part] = part_sum / nsmpl;
        }

        /* The average of the larger partitions is calculated by merging the averages of the smaller partitions */
        for (i = (int32_t)(max_porder - 1); i >= 0; i--) {
            for (part = 0; part < (1U << i); part++) {
                coder->part_mean[i][part] = (coder->part_mean[i + 1][2 * part] + coder->part_mean[i + 1][2 * part + 1]) / 2.0;
            }
        }
    }

    /* Switch signs based on overall average */
    if (max_uval == 0) {
        tmp_code_type = SRLACODER_CODE_TYPE_ALLZERO;
    } else if (coder->part_mean[0][0] < 2) {
        tmp_code_type = SRLACODER_CODE_TYPE_RICE;
    } else {
        tmp_code_type = SRLACODER_CODE_TYPE_RECURSIVE_RICE;
    }

    /* Calculate the code length for each division and search for the optimal division */
    min_bits = UINT32_MAX;
    best_porder = max_porder + 1;

    switch (tmp_code_type) {
    case SRLACODER_CODE_TYPE_ALLZERO:
        best_porder = 0;
        min_bits = 0;
        break;
    case SRLACODER_CODE_TYPE_RICE:
    {
        for (porder = 0; porder <= max_porder; porder++) {
            const uint32_t nsmpl = (num_samples >> porder);
            uint32_t k, prevk;
            uint32_t bits = SRLACODER_LOG2_MAX_NUM_PARTITIONS;
            for (part = 0; part < (1U << porder); part++) {
                SRLACoder_CalculateOptimalRiceParameter(coder->part_mean[porder][part], &k, NULL);
                for (smpl = 0; smpl < nsmpl; smpl++) {
                    bits += Rice_GetCodeLength(k, coder->uval_buffer[part * nsmpl + smpl]);
                }
                if (part == 0) {
                    bits += SRLACODER_RICE_PARAMETER_BITS;
                } else {
                    const int32_t diff = (int32_t)k - (int32_t)prevk;
                    const uint32_t udiff = SRLAUTILITY_SINT32_TO_UINT32(diff);
                    bits += udiff + 1;
                }
                prevk = k;
                /* If the minimum value is exceeded along the way, stop */
                if (bits >= min_bits) {
                    break;
                }
            }
            if (bits < min_bits) {
                min_bits = bits;
                best_porder = porder;
            }
        }
    }
    break;
    case SRLACODER_CODE_TYPE_RECURSIVE_RICE:
    {
        for (porder = 0; porder <= max_porder; porder++) {
            const uint32_t nsmpl = (num_samples >> porder);
            uint32_t k1, k2, prevk2;
            uint32_t bits = SRLACODER_LOG2_MAX_NUM_PARTITIONS;
            for (part = 0; part < (1U << porder); part++) {
                SRLACoder_CalculateOptimalRecursiveRiceParameter(coder->part_mean[porder][part], &k1, &k2, NULL);
                bits += RecursiveRice_ComputeCodeLength(&coder->uval_buffer[part * nsmpl], nsmpl, k1, k2);
                if (part == 0) {
                    bits += SRLACODER_RICE_PARAMETER_BITS;
                } else {
                    const int32_t diff = (int32_t)k2 - (int32_t)prevk2;
                    const uint32_t udiff = SRLAUTILITY_SINT32_TO_UINT32(diff);
                    bits += udiff + 1;
                }
                prevk2 = k2;
                /* If the minimum value is exceeded along the way, stop */
                if (bits >= min_bits) {
                    break;
                }
            }
            if (bits < min_bits) {
                min_bits = bits;
                best_porder = porder;
            }
        }
    }
    break;
    default:
        SRLA_ASSERT(0);
    }

    SRLA_ASSERT(best_porder != (max_porder + 1));

    /* Add signed type fields */
    min_bits += 2;

    /* Record the result */
    (*code_type) = tmp_code_type;
    (*best_partition_order) = best_porder;
    (*best_code_length) = min_bits;
}

/* Recursive Golomb--Rice encoding of data arrays */
static void SRLACoder_EncodeRecursiveRice(
    struct BitStream *stream, uint32_t *data, uint32_t num_samples, const uint32_t k1, const uint32_t k2)
{
    SRLA_ASSERT(stream != NULL);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(k1 == (k2 + 1));

    /* Define the decoding process for specific k2 parameters */
#define DEFINE_ENCODING_PROCEDURE_CASE(k2_param, stream, data, num_samples)\
    case (k2_param):\
        {\
            while ((num_samples)--) {\
                const uint32_t uval__ = (*((data)++));\
                RecursiveRice_PutCode((stream), ((k2_param)+1), (k2_param), uval__);\
            }\
        }\
        break;

    /* Case distinction by parameters */
    switch (k2) {
        DEFINE_ENCODING_PROCEDURE_CASE(  0, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  1, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  2, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  3, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  4, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  5, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  6, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  7, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  8, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE(  9, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 10, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 11, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 12, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 13, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 14, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 15, stream, data, num_samples);
        DEFINE_ENCODING_PROCEDURE_CASE( 16, stream, data, num_samples);
    default:
        while (num_samples--) {
            uint32_t uval = *(data++);
            RecursiveRice_PutCode(stream, k2 + 1, k2, uval);
        }
    }
}

/* Encoding a signed integer array */
static void SRLACoder_EncodePartitionedRecursiveRice(struct SRLACoder *coder, struct BitStream *stream, const int32_t *data, uint32_t num_samples)
{
    uint32_t part, best_porder, smpl, min_bits;
    SRLACoderCodeType code_type;

    /* Search for optimal code and partition */
    SRLACoder_SearchBestCodeAndPartition(
        coder, data, num_samples, &code_type, &best_porder, &min_bits);

    /* Encode using optimal partitioning */
    {
        SRLA_ASSERT(code_type != SRLACODER_CODE_TYPE_INVALID);
        BitWriter_PutBits(stream, code_type, 2);

        switch (code_type) {
        case SRLACODER_CODE_TYPE_ALLZERO:
            break;
        case SRLACODER_CODE_TYPE_RICE:
        {
            uint32_t k, prevk;
            const uint32_t nsmpl = num_samples >> best_porder;

            BitWriter_PutBits(stream, best_porder, SRLACODER_LOG2_MAX_NUM_PARTITIONS);
            for (part = 0; part < (1U << best_porder); part++) {
                SRLACoder_CalculateOptimalRiceParameter(coder->part_mean[best_porder][part], &k, NULL);
                if (part == 0) {
                    BitWriter_PutBits(stream, k, SRLACODER_RICE_PARAMETER_BITS);
                } else {
                    const int32_t diff = (int32_t)k - (int32_t)prevk;
                    const uint32_t udiff = SRLAUTILITY_SINT32_TO_UINT32(diff);
                    BitWriter_PutZeroRun(stream, udiff);
                }
                prevk = k;
                for (smpl = 0; smpl < nsmpl; smpl++) {
                    Rice_PutCode(stream, k, coder->uval_buffer[part * nsmpl + smpl]);
                }
            }
        }
        break;
        case SRLACODER_CODE_TYPE_RECURSIVE_RICE:
        {
            uint32_t k1, k2, prevk2;
            const uint32_t nsmpl = num_samples >> best_porder;

            BitWriter_PutBits(stream, best_porder, SRLACODER_LOG2_MAX_NUM_PARTITIONS);
            for (part = 0; part < (1U << best_porder); part++) {
                SRLACoder_CalculateOptimalRecursiveRiceParameter(coder->part_mean[best_porder][part], &k1, &k2, NULL);
                if (part == 0) {
                    BitWriter_PutBits(stream, k2, SRLACODER_RICE_PARAMETER_BITS);
                } else {
                    const int32_t diff = (int32_t)k2 - (int32_t)prevk2;
                    const uint32_t udiff = SRLAUTILITY_SINT32_TO_UINT32(diff);
                    BitWriter_PutZeroRun(stream, udiff);
                }
                prevk2 = k2;
                SRLACoder_EncodeRecursiveRice(stream, &coder->uval_buffer[part * nsmpl], nsmpl, k1, k2);
            }
        }
        break;
        default:
            SRLA_ASSERT(0);
        }
    }
}

/* Recursive Golomb--Rice decoding of data arrays */
static void SRLACoder_DecodeRecursiveRice(
    struct BitStream *stream, int32_t *data, uint32_t num_samples, const uint32_t k1, const uint32_t k2)
{
    SRLA_ASSERT(stream != NULL);
    SRLA_ASSERT(data != NULL);
    SRLA_ASSERT(k1 == (k2 + 1));

    /* Define the decoding process for specific k2 parameters */
#define DEFINE_DECODING_PROCEDURE_CASE(k2_param, stream, data, num_samples)\
    case (k2_param):\
        {\
            uint32_t uval__;\
            while ((num_samples)--) {\
                RecursiveRice_GetCode((stream), ((k2_param)+1), (k2_param), &uval__);\
                (*((data)++)) = SRLAUTILITY_UINT32_TO_SINT32(uval__);\
            }\
        }\
        break;

    /* Case distinction by parameters */
    switch (k2) {
        DEFINE_DECODING_PROCEDURE_CASE(  0, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  1, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  2, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  3, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  4, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  5, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  6, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  7, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  8, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE(  9, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 10, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 11, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 12, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 13, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 14, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 15, stream, data, num_samples);
        DEFINE_DECODING_PROCEDURE_CASE( 16, stream, data, num_samples);
    default:
        while (num_samples--) {
            uint32_t uval;
            RecursiveRice_GetCode(stream, k2 + 1, k2, &uval);
            *(data++) = SRLAUTILITY_UINT32_TO_SINT32(uval);
        }
    }
}

/* Decode signed integer array */
static void SRLACoder_DecodePartitionedRecursiveRice(struct BitStream *stream, int32_t *data, uint32_t num_samples)
{
    uint32_t smpl, part, nsmpl, best_porder;

    SRLACoderCodeType code_type = SRLACODER_CODE_TYPE_INVALID;

    BitReader_GetBits(stream, (uint32_t *)&code_type, 2);

    switch (code_type) {
    case SRLACODER_CODE_TYPE_ALLZERO:
        memset(data, 0, sizeof(int32_t) * num_samples);
        break;
    case SRLACODER_CODE_TYPE_RICE:
    {
        uint32_t k;
        BitReader_GetBits(stream, &best_porder, SRLACODER_LOG2_MAX_NUM_PARTITIONS);
        nsmpl = num_samples >> best_porder;
        for (part = 0; part < (1U << best_porder); part++) {
            if (part == 0) {
                BitReader_GetBits(stream, &k, SRLACODER_RICE_PARAMETER_BITS);
            } else {
                uint32_t udiff;
                BitReader_GetZeroRunLength(stream, &udiff);
                k = (uint32_t)((int32_t)k + SRLAUTILITY_UINT32_TO_SINT32(udiff));
            }
            for (smpl = 0; smpl < nsmpl; smpl++) {
                const uint32_t uval = Rice_GetCode(stream, k);
                data[part * nsmpl + smpl] = SRLAUTILITY_UINT32_TO_SINT32(uval);
            }
        }
    }
    break;
    case SRLACODER_CODE_TYPE_RECURSIVE_RICE:
    {
        uint32_t k2;
        BitReader_GetBits(stream, &best_porder, SRLACODER_LOG2_MAX_NUM_PARTITIONS);
        nsmpl = num_samples >> best_porder;
        for (part = 0; part < (1U << best_porder); part++) {
            if (part == 0) {
                BitReader_GetBits(stream, &k2, SRLACODER_RICE_PARAMETER_BITS);
            } else {
                uint32_t udiff;
                BitReader_GetZeroRunLength(stream, &udiff);
                k2 = (uint32_t)((int32_t)k2 + SRLAUTILITY_UINT32_TO_SINT32(udiff));
            }
            SRLACoder_DecodeRecursiveRice(stream, &data[part * nsmpl], nsmpl, k2 + 1, k2);
        }
    }
    break;
    default:
        SRLA_ASSERT(0);
    }
}

/* Code length calculation */
uint32_t SRLACoder_ComputeCodeLength(struct SRLACoder *coder, const int32_t *data, uint32_t num_samples)
{
    uint32_t dummy_best_porder, min_bits;
    SRLACoderCodeType dummy_code_type;

    SRLA_ASSERT((data != NULL) && (coder != NULL));
    SRLA_ASSERT(num_samples != 0);

    /* Search for optimal code and partition */
    SRLACoder_SearchBestCodeAndPartition(
        coder, data, num_samples, &dummy_code_type, &dummy_best_porder, &min_bits);

    return min_bits;
}

/* Encoding a signed integer array */
void SRLACoder_Encode(struct SRLACoder *coder, struct BitStream *stream, const int32_t *data, uint32_t num_samples)
{
    SRLA_ASSERT((stream != NULL) && (data != NULL) && (coder != NULL));
    SRLA_ASSERT(num_samples != 0);

    SRLACoder_EncodePartitionedRecursiveRice(coder, stream, data, num_samples);
}

/* Decode signed integer array */
void SRLACoder_Decode(struct BitStream *stream, int32_t *data, uint32_t num_samples)
{
    SRLA_ASSERT((stream != NULL) && (data != NULL));
    SRLA_ASSERT(num_samples != 0);

    SRLACoder_DecodePartitionedRecursiveRice(stream, data, num_samples);
}
