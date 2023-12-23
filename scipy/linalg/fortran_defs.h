/*
 * Handle different Fortran conventions.
 */

#include <stdlib.h>

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

#ifdef ACCELERATE_NEW_LAPACK
    #if __MAC_OS_X_VERSION_MAX_ALLOWED < 130300
        #ifdef HAVE_BLAS_ILP64
            #error "Accelerate ILP64 support is only available with macOS 13.3 SDK or later"
        #endif
    #else
        #define NO_APPEND_FORTRAN
        #ifdef HAVE_BLAS_ILP64
            #define BLAS_SYMBOL_SUFFIX $NEWLAPACK$ILP64
        #else
            #define BLAS_SYMBOL_SUFFIX $NEWLAPACK
        #endif
    #endif
#endif

#if defined(NO_APPEND_FORTRAN) || defined(ACCELERATE_NEW_LAPACK)
#define BLAS_FORTRAN_SUFFIX
#else
#define BLAS_FORTRAN_SUFFIX _
#endif

#ifndef BLAS_SYMBOL_PREFIX
#define BLAS_SYMBOL_PREFIX
#endif

#ifndef BLAS_SYMBOL_SUFFIX
#define BLAS_SYMBOL_SUFFIX
#endif

#define BLAS_FUNC_CONCAT(name,prefix,suffix,suffix2) prefix ## name ## suffix ## suffix2
#define BLAS_FUNC_EXPAND(name,prefix,suffix,suffix2) BLAS_FUNC_CONCAT(name,prefix,suffix,suffix2)

/*
 * Use either the OpenBLAS scheme with the `64_` suffix behind the Fortran
 * compiler symbol mangling, or the MKL scheme (and upcoming
 * reference-lapack#666) which does it the other way around and uses `_64`.
 */
#ifdef OPENBLAS_ILP64_NAMING_SCHEME
#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,BLAS_FORTRAN_SUFFIX,BLAS_SYMBOL_SUFFIX)
#else
#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,BLAS_SYMBOL_SUFFIX,BLAS_FORTRAN_SUFFIX)
#endif
/*
 * Note that CBLAS doesn't include Fortran compiler symbol mangling, so ends up
 * being the same in both schemes
 */
#define CBLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,,BLAS_SYMBOL_SUFFIX)

#ifdef HAVE_BLAS_ILP64
#define F_INT npy_int64
#else
#define F_INT int
#endif
