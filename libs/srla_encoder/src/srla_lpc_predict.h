#ifndef SRLA_LPCPREDICTOR_H_INCLUDED
#define SRLA_LPCPREDICTOR_H_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Prediction/error output by LPC coefficients */
void SRLALPC_Predict(
    const int32_t *data, uint32_t num_samples,
    const int32_t *coef, uint32_t coef_order, int32_t *residual, uint32_t coef_rshift);

#ifdef __cplusplus
}
#endif

#endif /* SRLA_LPCPREDICTOR_H_INCLUDED */
