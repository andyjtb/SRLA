/*
/*!
* @file fft.h
* @brief FFT (Fast Fourier Transform) library
*/
*/
#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

/*! @brief Access the real part of the i-th complex number */
#define FFTCOMPLEX_REAL(flt_array, i) ((flt_array)[((i) << 1)])
/*! @brief Access the imaginary part of the i-th complex number */
#define FFTCOMPLEX_IMAG(flt_array, i) ((flt_array)[((i) << 1) + 1])

#ifdef __cplusplus
extern "C" {
#endif

/*
/*!
* @brief FFT (fast Fourier transform)
* @param[in] n FFT points
* @param[in] flag -1: FFT, 1: IFFT
* @param[in,out] x Series to be Fourier transformed (input/output 2n size required, real part in even number, imaginary part in odd number)
* @param[in,out] y Working array (same size as x)
* @note No normalization is performed
*/
*/
void FFT_FloatFFT(int n, int flag, double *x, double *y);

/*
/*!
* @brief FFT (fast Fourier transform) of real array
* @param[in] n FFT points
* @param[in] flag -1: FFT, 1: IFFT
* @param[in,out] x Series to be Fourier transformed (input/output n size required, for FFT, x[0] contains the real part of the DC component, and x[1] contains the imaginary part of the highest frequency component)
* @param[in,out] y Working array (same size as x)
* @note No normalization is performed. The normalization constant is 2/n
*/
*/
void FFT_RealFFT(int n, int flag, double *x, double *y);

#ifdef __cplusplus
}
#endif

#endif /* FFT_H_INCLUDED */
