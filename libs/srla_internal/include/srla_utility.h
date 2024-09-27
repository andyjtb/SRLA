#ifndef SRLAUTILITY_H_INCLUDED
#define SRLAUTILITY_H_INCLUDED

#include "srla_stdint.h"
#include <stddef.h>

/* Avoid unused argument warning */
#define SRLAUTILITY_UNUSED_ARGUMENT(arg)  ((void)(arg))
/* Arithmetic right shift */
#if ((((int32_t)-1) >> 1) == ((int32_t)-1))
/* In an environment where arithmetic right shift is valid, it is simply shifted right */
#define SRLAUTILITY_SHIFT_RIGHT_ARITHMETIC(sint32, rshift) ((sint32) >> (rshift))
#else
/* In environments where arithmetic right shift is disabled, define it yourself. Quote from Hacker's Fun */
/* Note) Valid range: 0 <= rshift <= 32 */
#define SRLAUTILITY_SHIFT_RIGHT_ARITHMETIC(sint32, rshift) ((((uint64_t)(sint32) + 0x80000000UL) >> (rshift)) - (0x80000000UL >> (rshift)))
#endif
/* Sign function. Cited from Hacker's Fun. Supplementary Note: Returns 0 when val==0. */
#define SRLAUTILITY_SIGN(val)  (((val) > 0) - ((val) < 0))
/* Round up to the next multiple of n */
#define SRLAUTILITY_ROUNDUP(val, n) ((((val) + ((n) - 1)) / (n)) * (n))
/* Get the maximum value */
#define SRLAUTILITY_MAX(a,b) (((a) > (b)) ? (a) : (b))
/* Get the minimum value */
#define SRLAUTILITY_MIN(a,b) (((a) < (b)) ? (a) : (b))
/* Limit to minimum value or more and minimum value or less */
#define SRLAUTILITY_INNER_VALUE(val, min, max) (SRLAUTILITY_MIN((max), SRLAUTILITY_MAX((min), (val))))
/* Is it a power of 2? */
#define SRLAUTILITY_IS_POWERED_OF_2(val) (!((val) & ((val) - 1)))
/* Uniquely convert a signed 32-bit number to an unsigned 32-bit number */
#define SRLAUTILITY_SINT32_TO_UINT32(sint) ((uint32_t)((-((sint) < 0)) ^ ((sint) << 1)))
/* Uniquely convert an unsigned 32-bit number to a signed 32-bit number */
#define SRLAUTILITY_UINT32_TO_SINT32(uint) ((int32_t)((uint) >> 1) ^ -(int32_t)((uint) & 1))
/* Get absolute value */
#define SRLAUTILITY_ABS(val)               (((val) > 0) ? (val) : -(val))

/* Calculate NLZ (number of bits from the most significant bit to 1) */
#if defined(__GNUC__)
/* Use built-in functions */
#define SRLAUTILITY_NLZ(x) (((x) > 0) ? (uint32_t)__builtin_clz(x) : 32U)
#elif defined(_MSC_VER)
#include <intrin.h>
/* Use built-in functions */
__inline uint32_t SRLAUTILITY_NLZ(uint32_t x)
{
    return __lzcnt(x);
}
#else
/* Use software implementation */
#define SRLAUTILITY_NLZ(x) SRLAUtility_NLZSoft(x)
#endif

/* Calculate ceil(log2(val)) */
#define SRLAUTILITY_LOG2CEIL(x) (32U - SRLAUTILITY_NLZ((uint32_t)((x) - 1U)))
/* Calculate floor(log2(val)) */
#define SRLAUTILITY_LOG2FLOOR(x) (31U - SRLAUTILITY_NLZ(x))

/* Round up to a power of 2 (1,2,4,8,16,...) */
#if defined(__GNUC__) || defined(_MSC_VER)
/* Use built-in functions */
#define SRLAUTILITY_ROUNDUP2POWERED(x) (1U << SRLAUTILITY_LOG2CEIL(x))
#else
/* Use software implementation */
#define SRLAUTILITY_ROUNDUP2POWERED(x) SRLAUtility_RoundUp2PoweredSoft(x)
#endif

/* Calculate the area work size of a two-dimensional array */
#define SRLA_CALCULATE_2DIMARRAY_WORKSIZE(type, size1, size2)\
    ((size1) * ((int32_t)sizeof(type *) + SRLA_MEMORY_ALIGNMENT\
        + (size2) * (int32_t)sizeof(type) + SRLA_MEMORY_ALIGNMENT))

/* Allocate space for 2-dimensional array */
#define SRLA_ALLOCATE_2DIMARRAY(ptr, work_ptr, type, size1, size2)\
    do {\
        uint32_t i;\
        (work_ptr) = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);\
        (ptr) = (type **)work_ptr;\
        (work_ptr) += sizeof(type *) * (size1);\
        for (i = 0; i < (size1); i++) {\
            (work_ptr) = (uint8_t *)SRLAUTILITY_ROUNDUP((uintptr_t)work_ptr, SRLA_MEMORY_ALIGNMENT);\
            (ptr)[i] = (type *)work_ptr;\
            (work_ptr) += sizeof(type) * (size2);\
        }\
    } while (0)

/* Pre-emphasis/de-emphasis filter */
struct SRLAPreemphasisFilter {
    int32_t prev;
    int32_t coef;
};

#ifdef __cplusplus
extern "C" {
#endif

/* round function (not provided in C89) */
double SRLAUtility_Round(double d);

/* log2 function (not defined in C89) */
double SRLAUtility_Log2(double d);

/* Fletcher checksum calculation */
uint16_t SRLAUtility_CalculateFletcher16CheckSum(const uint8_t *data, size_t data_size);

/* Calculate NLZ (number of bits from the most significant bit to 1) */
uint32_t SRLAUtility_NLZSoft(uint32_t val);

/* Round up to the next power of 2 */
uint32_t SRLAUtility_RoundUp2PoweredSoft(uint32_t val);

/* LR -> MS (in-place) */
void SRLAUtility_LRtoMSConversion(int32_t **buffer, uint32_t num_samples);

/* MS -> LR (in-place) */
void SRLAUtility_MStoLRConversion(int32_t **buffer, uint32_t num_samples);

/* LR -> LS (in-place) */
void SRLAUtility_LRtoLSConversion(int32_t **buffer, uint32_t num_samples);

/* LS -> LR (in-place) */
void SRLAUtility_LStoLRConversion(int32_t **buffer, uint32_t num_samples);

/* LR -> SR (in-place) */
void SRLAUtility_LRtoSRConversion(int32_t **buffer, uint32_t num_samples);

/* SR -> LR (in-place) */
void SRLAUtility_SRtoLRConversion(int32_t **buffer, uint32_t num_samples);

/* Pre-emphasis filter initialization */
void SRLAPreemphasisFilter_Initialize(struct SRLAPreemphasisFilter *preem);

/* Calculate the coefficient of multi-stage pre-emphasis */
void SRLAPreemphasisFilter_CalculateMultiStageCoefficients(
    struct SRLAPreemphasisFilter *preem, uint32_t num_preem, const int32_t *buffer, uint32_t num_samples);

/* Pre-emphasis */
void SRLAPreemphasisFilter_Preemphasis(
    struct SRLAPreemphasisFilter *preem, int32_t *buffer, uint32_t num_samples);

/* Apply de-emphasis multiple times */
void SRLAPreemphasisFilter_MultiStageDeemphasis(
    struct SRLAPreemphasisFilter *preem, uint32_t num_preem, int32_t *buffer, uint32_t num_samples);

#ifdef __cplusplus
}
#endif

#endif /* SRLAUTILITY_H_INCLUDED */
