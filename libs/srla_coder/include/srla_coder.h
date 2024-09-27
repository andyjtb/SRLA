#ifndef SRLACODER_H_INCLUDED
#define SRLACODER_H_INCLUDED

#include <stdint.h>
#include "bit_stream.h"

/* encoding handle */
struct SRLACoder;

#ifdef __cplusplus
extern "C" {
#endif

/* Calculate the work size required to create the encoding handle */
int32_t SRLACoder_CalculateWorkSize(uint32_t max_num_samples);

/* Create an encoding handle */
struct SRLACoder* SRLACoder_Create(uint32_t max_num_samples, void *work, int32_t work_size);

/* Destroy the encoding handle */
void SRLACoder_Destroy(struct SRLACoder *coder);

/* Code length calculation */
uint32_t SRLACoder_ComputeCodeLength(struct SRLACoder *coder, const int32_t *data, uint32_t num_samples);

/* Encoding a signed integer array */
void SRLACoder_Encode(struct SRLACoder *coder, struct BitStream *stream, const int32_t *data, uint32_t num_samples);

/* Decode signed integer array */
void SRLACoder_Decode(struct BitStream *stream, int32_t *data, uint32_t num_samples);

#ifdef __cplusplus
}
#endif

#endif /* SRLACODER_H_INCLUDED */
