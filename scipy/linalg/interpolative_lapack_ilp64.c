#include "fortran_defs.h"
#include "numpy/arrayobject.h"

void BLAS_FUNC(dgesdd)(char *jobz, F_INT *m, F_INT *n, double *a, F_INT *lda, double *s, double *u, F_INT *ldu, double *vt, F_INT *ldvt, double *work, F_INT *lwork, F_INT *iwork, F_INT *info);
void BLAS_FUNC(zgesdd)(char *jobz, F_INT *m, F_INT *n, npy_complex128 *a, F_INT *lda, double *s, npy_complex128 *u, F_INT *ldu, npy_complex128 *vt, F_INT *ldvt, npy_complex128 *work, F_INT *lwork, double *rwork, F_INT *iwork, F_INT *info);

#ifdef HAVE_BLAS_ILP64
void F_FUNC(wdgesdd, WDGESDD)(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *iwork, int *info){
    # iwork has minimum size 8*min(m,n)
    int iwork_dim = 8 * *m;
    F_INT* iwork_temp = malloc(iwork_dim * sizeof(F_INT));
    if (iwork_temp == NULL) {
        fprintf(stderr, "Memory allocation for temporary 64-bit int array failed in wzgesdd\n");
        free(iwork_temp);
        exit(EXIT_FAILURE);
    }
    F_INT m_temp = *m;
    F_INT n_temp = *n;
    F_INT lda_temp = *lda;
    F_INT ldu_temp = *ldu;
    F_INT ldvt_temp = *ldvt;
    F_INT lwork_temp = *lwork;
    F_INT info_temp = *info;
    BLAS_FUNC(dgesdd)(jobz,&m_temp,&n_temp,a,&lda_temp,s,u,&ldu_temp,vt,&ldvt_temp,work,&lwork_temp,iwork_temp,&info_temp);
    # info and iwork are the only variables mutated by zgesdd
    *info = (int)info_temp;
    for (size_t i = 0; i < iwork_dim; ++i) {
        iwork[i] = (int)(iwork_temp[i]);
    }
    free(iwork_temp);
}

void F_FUNC(wzgesdd, WZGESDD)(char *jobz, int *m, int *n, npy_complex128 *a, int *lda, double *s, npy_complex128 *u, int *ldu, npy_complex128 *vt, int *ldvt, npy_complex128 *work, int *lwork, double *rwork, int *iwork, int *info){
    # iwork has minimum size 8*min(m,n)
    int iwork_dim = 8* *m;
    F_INT* iwork_temp = malloc(iwork_dim * sizeof(F_INT));
    if (iwork_temp == NULL) {
        fprintf(stderr, "Memory allocation for temporary 64-bit int array failed in wzgesdd\n");
        free(iwork_temp);
        exit(EXIT_FAILURE);
    }
    F_INT m_temp = *m;
    F_INT n_temp = *n;
    F_INT lda_temp = *lda;
    F_INT ldu_temp = *ldu;
    F_INT ldvt_temp = *ldvt;
    F_INT lwork_temp = *lwork;
    F_INT info_temp = *info;
    BLAS_FUNC(zgesdd)(jobz,&m_temp,&n_temp,a,&lda_temp,s,u,&ldu_temp,vt,&ldvt_temp,work,&lwork_temp,rwork,iwork_temp,&info_temp);
    # info and iwork are the only variables mutated by zgesdd
    *info = (int)info_temp;
    for (size_t i = 0; i < iwork_dim; ++i) {
        iwork[i] = (int)(iwork_temp[i]);
    }
    free(iwork_temp);
}

#else
void F_FUNC(wdgesdd, WDGESDD)(char *jobz, F_INT *m, F_INT *n, double *a, F_INT *lda, double *s, double *u, F_INT *ldu, double *vt, F_INT *ldvt, double *work, F_INT *lwork, F_INT *iwork, F_INT *info){
    BLAS_FUNC(dgesdd)(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info);
}
void F_FUNC(wzgesdd, WZGESDD)(char *jobz, F_INT *m, F_INT *n, npy_complex128 *a, F_INT *lda, double *s, npy_complex128 *u, F_INT *ldu, npy_complex128 *vt, F_INT *ldvt, npy_complex128 *work, F_INT *lwork, double *rwork, F_INT *iwork, F_INT *info){
    BLAS_FUNC(zgesdd)(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,iwork,info);
}

#endif
