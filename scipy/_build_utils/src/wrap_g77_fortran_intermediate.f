!      SUBROUTINE wrappers for BLAS/LAPACK FUNCTIONs with
!      complex return values. Fixes x86 segfaults due to
!      incompatibility between struct and C99 complex types.
       SUBROUTINE CDOTCWRP( RET, N, CX, INCX, CY, INCY )
       INTEGER INCX, INCY, N
       COMPLEX CX(*), CY(*), RET
       EXTERNAL WCDOTC
       COMPLEX WCDOTC
       RET = WCDOTC( N, CX, INCX, CY, INCY )
       END SUBROUTINE

       SUBROUTINE CDOTUWRP( RET, N, CX, INCX, CY, INCY )
       INTEGER INCX, INCY, N
       COMPLEX CX(*), CY(*), RET
       EXTERNAL WCDOTU
       COMPLEX WCDOTU
       RET = WCDOTU( N, CX, INCX, CY, INCY )
       END SUBROUTINE

       SUBROUTINE ZDOTCWRP( RET, N, CX, INCX, CY, INCY )
       INTEGER INCX, INCY, N
       DOUBLE COMPLEX CX(*), CY(*), RET
       EXTERNAL WZDOTC
       DOUBLE COMPLEX WZDOTC
       RET = WZDOTC( N, CX, INCX, CY, INCY )
       END SUBROUTINE

       SUBROUTINE ZDOTUWRP( RET, N, CX, INCX, CY, INCY )
       INTEGER INCX, INCY, N
       DOUBLE COMPLEX CX(*), CY(*), RET
       EXTERNAL WZDOTU
       DOUBLE COMPLEX WZDOTU
       RET = WZDOTU( N, CX, INCX, CY, INCY )
       END SUBROUTINE

       SUBROUTINE CLADIVWRP( RET, X, Y )
       COMPLEX            X, Y, RET
       EXTERNAL           WCLADIV
       COMPLEX            WCLADIV
       RET = WCLADIV( X, Y)
       END SUBROUTINE

       SUBROUTINE ZLADIVWRP( RET, X, Y )
       DOUBLE COMPLEX     X, Y, RET
       EXTERNAL           WZLADIV
       DOUBLE COMPLEX     WZLADIV
       RET = WZLADIV( X, Y)
       END SUBROUTINE
