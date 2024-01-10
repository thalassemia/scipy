#include "fortran_defs.h"
#include "npy_cblas.h"

#ifdef F_SRC
#include <complex.h>
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
typedef _Dcomplex double_complex;
typedef _Fcomplex float_complex;
#else /* !defined(_MSC_VER) || defined(__INTEL_COMPILER) */
typedef double _Complex double_complex;
typedef float _Complex float_complex;
#endif
#else
typedef npy_complex128 double_complex;
typedef npy_complex64 float_complex;
#endif

float_complex F_FUNC(wcdotc, WCDOTC)(CBLAS_INT *n, float_complex *cx, CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    float_complex ret;
    CBLAS_FUNC(cblas_cdotc_sub)(*n, cx, *incx, cy, *incy,&ret);
    return ret;
}

double_complex F_FUNC(wzdotc, WZDOTC)(CBLAS_INT *n, double_complex *zx, CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    double_complex ret;
    CBLAS_FUNC(cblas_zdotc_sub)(*n, zx, *incx, zy, *incy,&ret);
    return ret;
}

float_complex F_FUNC(wcdotu, WCDOTU)(CBLAS_INT *n, float_complex *cx, CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    float_complex ret;
    CBLAS_FUNC(cblas_cdotu_sub)(*n, cx, *incx, cy, *incy,&ret);
    return ret;
}

double_complex F_FUNC(wzdotu, WZDOTU)(CBLAS_INT *n, double_complex *zx, CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    double_complex ret;
    CBLAS_FUNC(cblas_zdotu_sub)(*n, zx, *incx, zy, *incy,&ret);
    return ret;
}

void BLAS_FUNC(sladiv)(float *xr, float *xi, float *yr, float *yi, float *retr, float * reti);
void F_FUNC(wcladiv, WCLADIV)(void *ret, void *x, void *y){
    BLAS_FUNC(sladiv)((float*)(x), (float*)(x)+1, \
        (float*)(y), (float*)(y)+1, \
        (float*)(ret), (float*)(ret)+1);
}

void BLAS_FUNC(dladiv)(double *xr, double *xi, double *yr, double *yi, double *retr, double * reti);
void F_FUNC(wzladiv, WZLADIV)(void *ret, void *x, void *y){
    BLAS_FUNC(dladiv)((double*)(x), (double*)(x)+1, \
        (double*)(y), (double*)(y)+1, \
        (double*)(ret), (double*)(ret)+1);
}
