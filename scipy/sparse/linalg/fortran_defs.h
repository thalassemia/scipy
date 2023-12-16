/*
 * Handle different Fortran conventions.
 */

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
    #define BLAS_SYMBOL_SUFFIX $NEWLAPACK
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

#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name, BLAS_SYMBOL_PREFIX, BLAS_SYMBOL_SUFFIX, BLAS_FORTRAN_SUFFIX)
