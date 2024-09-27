#include "fft.h"

#include <string.h>
#include <math.h>
#include <assert.h>

/* Define inline keywords */
#if defined(_MSC_VER)
#define FFT_INLINE inline
#elif defined(__GNUC__)
#define FFT_INLINE __inline__
#else
#define FFT_INLINE
#endif

/* Pi */
#define FFT_PI 3.14159265358979323846

/* Complex number type */
typedef struct FFTComplex {
    double real; /* real part */
    double imag; /* imaginary part */
} FFTComplex;

/*
/* FFT normalization is not performed
* n sequence length
* flag -1: FFT, 1: IFFT
* x sequence to be Fourier transformed (input/output)
* y work array (same size as x)
*/
*/
static void FFT_ComplexFFT(int n, int flag, FFTComplex *x, FFTComplex *y);

/*
/* Check the size of complex numbers. Double arrays are treated as complex numbers for calculations.
* If padding is added to the structure, the size will not match.
* If it does not match, pack the structure with #pragma. */
*/
extern char FFT_checksize[(sizeof(FFTComplex) == (sizeof(double) * 2)) ? 1 : -1];

/* Complex addition */
static FFT_INLINE FFTComplex FFTComplex_Add(FFTComplex a, FFTComplex b)
{
    FFTComplex ret;
    ret.real = a.real + b.real;
    ret.imag = a.imag + b.imag;
    return ret;
}

/* Complex subtraction */
static FFT_INLINE FFTComplex FFTComplex_Sub(FFTComplex a, FFTComplex b)
{
    FFTComplex ret;
    ret.real = a.real - b.real;
    ret.imag = a.imag - b.imag;
    return ret;
}

/* Complex multiplication */
static FFT_INLINE FFTComplex FFTComplex_Mul(FFTComplex a, FFTComplex b)
{
    FFTComplex ret;
    ret.real = a.real * b.real - a.imag * b.imag;
    ret.imag = a.real * b.imag + a.imag * b.real;
    return ret;
}

/*
/* FFT normalization is not performed
* n sequence length
* flag -1: FFT, 1: IFFT
* x sequence to be Fourier transformed (input/output)
* y work array (same size as x)
*/
*/
static void FFT_ComplexFFT(int n, const int flag, FFTComplex *x, FFTComplex *y)
{
    FFTComplex *tmp, *src = x;
    int p, q;
    int s = 1; /* stride */

    /* radix 4 Stockham FFT */
    while (n > 2) {
        const int n1 = (n >> 2);
        const int n2 = (n >> 1);
        const int n3 = n1 + n2;
        const double theta0 = 2.0 * FFT_PI / n;
        FFTComplex j, wdelta, w1p;
        j.real = 0.0; j.imag = -flag;
        wdelta.real = cos(theta0); wdelta.imag = flag * sin(theta0);
        w1p.real = 1.0; w1p.imag = 0.0;
        for (p = 0; p < n1; p++) {
            /*
/* More precise, but has sin and cos function calls. * const FFTComplex w1p = { cos(p * theta0), flag * sin(p * theta0) }; */
*/
            const FFTComplex w2p = FFTComplex_Mul(w1p, w1p);
            const FFTComplex w3p = FFTComplex_Mul(w1p, w2p);
            for (q = 0; q < s; q++) {
                const FFTComplex    a = x[q + s * (p +  0)];
                const FFTComplex    b = x[q + s * (p + n1)];
                const FFTComplex    c = x[q + s * (p + n2)];
                const FFTComplex    d = x[q + s * (p + n3)];
                const FFTComplex  apc = FFTComplex_Add(a, c);
                const FFTComplex  amc = FFTComplex_Sub(a, c);
                const FFTComplex  bpd = FFTComplex_Add(b, d);
                const FFTComplex jbmd = FFTComplex_Mul(j, FFTComplex_Sub(b, d));
                y[q + s * ((p << 2) + 0)] = FFTComplex_Add(apc, bpd);
                y[q + s * ((p << 2) + 1)] = FFTComplex_Mul(w1p, FFTComplex_Sub(amc, jbmd));
                y[q + s * ((p << 2) + 2)] = FFTComplex_Mul(w2p, FFTComplex_Sub(apc,  bpd));
                y[q + s * ((p << 2) + 3)] = FFTComplex_Mul(w3p, FFTComplex_Add(amc, jbmd));
            }
            /* Advance rotation factor */
            w1p = FFTComplex_Mul(w1p, wdelta);
        }
        n >>= 2;
        s <<= 2;
        tmp = x; x = y; y = tmp;
    }

    if (n == 2) {
        for (q = 0; q < s; q++) {
            const FFTComplex a = x[q + 0];
            const FFTComplex b = x[q + s];
            y[q + 0] = FFTComplex_Add(a, b);
            y[q + s] = FFTComplex_Sub(a, b);
        }
        s <<= 1;
        tmp = x; x = y; y = tmp;
    }

    if (src != x) {
        memcpy(y, x, sizeof(FFTComplex) * (size_t)s);
    }
}

/*
/* FFT normalization is not performed
* n sequence length
* flag -1: FFT, 1: IFFT
* x sequence to be Fourier transformed (input/output size 2n required, real part in even numbers, imaginary part in odd numbers)
* y work array (same size as x)
*/
*/
void FFT_FloatFFT(int n, const int flag, double *x, double *y)
{
    FFT_ComplexFFT(n, flag, (FFTComplex *)x, (FFTComplex *)y);
}

/*
/* FFT of real sequence. No normalization is performed. Normalization constant is 2/n
* n sequence length
* flag -1: FFT, 1: IFFT
* x sequence to be Fourier transformed (input/output size n required, in case of FFT, x[0] contains real part of DC component, x[1] contains imaginary part of highest frequency component)
* y work array (same size as x)
*/
*/
void FFT_RealFFT(int n, const int flag, double *x, double *y)
{
    int i;
    const double theta = flag * 2.0 * FFT_PI / n;
    const double wpi = sin(theta);
    const double wpr = cos(theta) - 1.0;
    const double c2 = flag * 0.5;
    double wr, wi, wtmp;

    /* In case of FFT, convert first */
    if (flag == -1) {
        FFT_FloatFFT(n >> 1, -1, x, y);
    }

    /* Rotation factor initialization */
    wr = 1.0 + wpr;
    wi = wpi;

    /* Using spectral symmetry */
    /* In the case of FFT, we put together the final results, and in the case of IFFT, we rearrange them back to their original state */
    for (i = 1; i <= (n >> 2); i++) {
        const int i1 = (i << 1);
        const int i2 = i1 + 1;
        const int i3 = n - i1;
        const int i4 = i3 + 1;
        const double h1r = 0.5 * (x[i1] + x[i3]);
        const double h1i = 0.5 * (x[i2] - x[i4]);
        const double h2r = -c2 * (x[i2] + x[i4]);
        const double h2i =  c2 * (x[i1] - x[i3]);
        x[i1] =  h1r + (wr * h2r) - (wi * h2i);
        x[i2] =  h1i + (wr * h2i) + (wi * h2r);
        x[i3] =  h1r - (wr * h2r) + (wi * h2i);
        x[i4] = -h1i + (wr * h2i) + (wi * h2r);
        /* Update twiddle factors */
        wtmp = wr;
        wr += wtmp * wpr - wi * wpi;
        wi += wi * wpr + wtmp * wpi;
    }

    /* DC component/highest frequency component */
    {
        const double h1r = x[0];
        if (flag == -1) {
            x[0] = h1r + x[1];
            x[1] = h1r - x[1];
        } else {
            x[0] = 0.5 * (h1r + x[1]);
            x[1] = 0.5 * (h1r - x[1]);
            FFT_FloatFFT(n >> 1, 1, x, y);
        }
    }
}
