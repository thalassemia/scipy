"""
Code generator script to make the Cython BLAS and LAPACK wrappers
from the files "cython_blas_signatures.txt" and
"cython_lapack_signatures.txt" which contain the signatures for
all the BLAS/LAPACK routines that should be included in the wrappers.
"""

import os
from stat import ST_MTIME
import argparse


BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Used to convert from types in signature files to Fortran types
fortran_types = {'int': 'integer',
                 'c': 'complex',
                 'd': 'double precision',
                 's': 'real',
                 'z': 'complex*16',
                 'char': 'character',
                 'bint': 'logical'}

# Used to convert from types in signature files to C types
c_types = {'int': 'int',
           'c': 'npy_complex64',
           'd': 'double',
           's': 'float',
           'z': 'npy_complex128',
           'char': 'char',
           'bint': 'int',
           'void': 'void',
           'cselect1': '_cselect1',
           'cselect2': '_cselect2',
           'dselect2': '_dselect2',
           'dselect3': '_dselect3',
           'sselect2': '_sselect2',
           'sselect3': '_sselect3',
           'zselect1': '_zselect1',
           'zselect2': '_zselect2'}

# Used to convert complex types in signature files to Numpy complex types
npy_types = {'c': 'npy_complex64', 'z': 'npy_complex128',
             'cselect1': '_cselect1', 'cselect2': '_cselect2',
             'dselect2': '_dselect2', 'dselect3': '_dselect3',
             'sselect2': '_sselect2', 'sselect3': '_sselect3',
             'zselect1': '_zselect1', 'zselect2': '_zselect2'}

# BLAS/LAPACK functions with complex return values (use 'wrp'-suffixed
# wrappers from G77 ABI wrapper)
wrapped_funcs = ['cdotc', 'cdotu', 'zdotc', 'zdotu', 'cladiv', 'zladiv']


def read_signatures(lines):
    """
    Read BLAS/LAPACK signatures and split into name, return type, C parameter
    list (e.g. `int* a, bint* b`), and argument names (e.g. `a, b`).
    """
    sigs = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        line = line[:-1].split('(')
        args = line[1]
        name_and_type = line[0].split(' ')
        ret_type = name_and_type[0]
        name = name_and_type[1]
        argtypes, argnames = zip(*[arg.split(' *') for arg in args.split(', ')])
        # Argname cannot be same as abbreviated return type
        if ret_type in argnames:
            argnames = [n if n != ret_type else n + '_' for n in argnames]
        # Argname should not be Python keyword
        argnames = [n if n not in ['lambda', 'in'] else n + '_' for n in argnames]
        sigs.append({
            'name': name,
            'return_type': ret_type,
            'argnames': argnames,
            'argtypes': argtypes
        })
    return sigs


def arg_casts(arg):
    """Cast from Cython to Numpy complex pointer type."""
    if arg in npy_types.values():
        return f'<{arg}*>'
    return ''


def pyx_decl(name, return_type, argnames, argtypes, header_name, accelerate):
    """Create Cython declaration for BLAS/LAPACK function."""
    args = ', '.join([' *'.join(arg) for arg in zip(argtypes, argnames)])
    argtypes = [npy_types.get(t, t) for t in argtypes]
    fort_args = ', '.join([' *'.join(arg) for arg in zip(argtypes, argnames)])
    argnames = ', '.join([arg_casts(t) + n for n, t in zip(argnames, argtypes)])

    fort_name = name
    fort_macro = 'BLAS_FUNC'
    # For functions with complex return type, use 'wrp'-suffixed wrappers
    if name in wrapped_funcs:
        # Cast from Cython to Numpy complex types
        npy_ret_type = npy_types.get(return_type, return_type)
        return f"""
cdef extern from "{header_name}":
    void _fortran_{name} "{name}wrp_"({npy_ret_type} *out, {fort_args}) nogil
cdef {return_type} {name}({args}) noexcept nogil:
    cdef {return_type} out
    _fortran_{name}(<{npy_ret_type}*>&out, {argnames})
    return out
"""
    # Only return if not void function
    return_kw = 'return ' if return_type != 'void' else ''
    # Not in new Accelerate (macOS 13.3+) so fallback to old
    if accelerate:
        if name in ['lsame', 'dcabs1']:
            fort_macro = ''
            fort_name = f'{name}_'
        elif name == 'xerbla_array':
            fort_macro = ''
            fort_name += '__'
    return f"""
cdef extern from "{header_name}":
    {return_type} _fortran_{name} "{fort_macro}({fort_name})"({fort_args}) nogil
cdef {return_type} {name}({args}) noexcept nogil:
    {return_kw}_fortran_{name}({argnames})
"""


blas_pyx_preamble = '''# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True

"""
BLAS Functions for Cython
=========================

Usable from Cython via::

    cimport scipy.linalg.cython_blas

These wrappers do not check for alignment of arrays.
Alignment should be checked before these wrappers are used.

Raw function pointers (Fortran-style pointer arguments):

- {}


"""

# Within SciPy, these wrappers can be used via relative or absolute cimport.
# Examples:
# from ..linalg cimport cython_blas
# from scipy.linalg cimport cython_blas
# cimport scipy.linalg.cython_blas as cython_blas
# cimport ..linalg.cython_blas as cython_blas

# Within SciPy, if BLAS functions are needed in C/C++/Fortran,
# these wrappers should not be used.
# The original libraries should be linked directly.

cdef extern from "fortran_defs.h":
    pass

from numpy cimport npy_complex64, npy_complex128

'''

lapack_pyx_preamble = '''"""
LAPACK functions for Cython
===========================

Usable from Cython via::

    cimport scipy.linalg.cython_lapack

This module provides Cython-level wrappers for all primary routines included
in LAPACK 3.4.0 except for ``zcgesv`` since its interface is not consistent
from LAPACK 3.4.0 to 3.6.0. It also provides some of the
fixed-api auxiliary routines.

These wrappers do not check for alignment of arrays.
Alignment should be checked before these wrappers are used.

Raw function pointers (Fortran-style pointer arguments):

- {}


"""

# Within SciPy, these wrappers can be used via relative or absolute cimport.
# Examples:
# from ..linalg cimport cython_lapack
# from scipy.linalg cimport cython_lapack
# cimport scipy.linalg.cython_lapack as cython_lapack
# cimport ..linalg.cython_lapack as cython_lapack

# Within SciPy, if LAPACK functions are needed in C/C++/Fortran,
# these wrappers should not be used.
# The original libraries should be linked directly.

cdef extern from "fortran_defs.h":
    pass

from numpy cimport npy_complex64, npy_complex128

cdef extern from "_lapack_subroutines.h":
    # Function pointer type declarations for
    # gees and gges families of functions.
    ctypedef bint _cselect1(npy_complex64*)
    ctypedef bint _cselect2(npy_complex64*, npy_complex64*)
    ctypedef bint _dselect2(d*, d*)
    ctypedef bint _dselect3(d*, d*, d*)
    ctypedef bint _sselect2(s*, s*)
    ctypedef bint _sselect3(s*, s*, s*)
    ctypedef bint _zselect1(npy_complex128*)
    ctypedef bint _zselect2(npy_complex128*, npy_complex128*)

'''

blas_py_wrappers = """

# Python-accessible wrappers for testing:

cdef inline bint _is_contiguous(double[:,:] a, int axis) noexcept nogil:
    return (a.strides[axis] == sizeof(a[0,0]) or a.shape[axis] == 1)

cpdef float complex _test_cdotc(float complex[:] cx, float complex[:] cy) noexcept nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return cdotc(&n, &cx[0], &incx, &cy[0], &incy)

cpdef float complex _test_cdotu(float complex[:] cx, float complex[:] cy) noexcept nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return cdotu(&n, &cx[0], &incx, &cy[0], &incy)

cpdef double _test_dasum(double[:] dx) noexcept nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return dasum(&n, &dx[0], &incx)

cpdef double _test_ddot(double[:] dx, double[:] dy) noexcept nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
        int incy = dy.strides[0] // sizeof(dy[0])
    return ddot(&n, &dx[0], &incx, &dy[0], &incy)

cpdef int _test_dgemm(double alpha, double[:,:] a, double[:,:] b, double beta,
                double[:,:] c) except -1 nogil:
    cdef:
        char *transa
        char *transb
        int m, n, k, lda, ldb, ldc
        double *a0=&a[0,0]
        double *b0=&b[0,0]
        double *c0=&c[0,0]
    # In the case that c is C contiguous, swap a and b and
    # swap whether or not each of them is transposed.
    # This can be done because a.dot(b) = b.T.dot(a.T).T.
    if _is_contiguous(c, 1):
        if _is_contiguous(a, 1):
            transb = 'n'
            ldb = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif _is_contiguous(a, 0):
            transb = 't'
            ldb = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if _is_contiguous(b, 1):
            transa = 'n'
            lda = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif _is_contiguous(b, 0):
            transa = 't'
            lda = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        k = b.shape[0]
        if k != a.shape[1]:
            with gil:
                raise ValueError("Shape mismatch in input arrays.")
        m = b.shape[1]
        n = a.shape[0]
        if n != c.shape[0] or m != c.shape[1]:
            with gil:
                raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[1,0]) - c0 if c.shape[0] > 1 else 1
        dgemm(transa, transb, &m, &n, &k, &alpha, b0, &lda, a0,
                   &ldb, &beta, c0, &ldc)
    elif _is_contiguous(c, 0):
        if _is_contiguous(a, 1):
            transa = 't'
            lda = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif _is_contiguous(a, 0):
            transa = 'n'
            lda = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if _is_contiguous(b, 1):
            transb = 't'
            ldb = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif _is_contiguous(b, 0):
            transb = 'n'
            ldb = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        m = a.shape[0]
        k = a.shape[1]
        if k != b.shape[0]:
            with gil:
                raise ValueError("Shape mismatch in input arrays.")
        n = b.shape[1]
        if m != c.shape[0] or n != c.shape[1]:
            with gil:
                raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[0,1]) - c0 if c.shape[1] > 1 else 1
        dgemm(transa, transb, &m, &n, &k, &alpha, a0, &lda, b0,
                   &ldb, &beta, c0, &ldc)
    else:
        with gil:
            raise ValueError("Input 'c' is neither C nor Fortran contiguous.")
    return 0

cpdef double _test_dnrm2(double[:] x) noexcept nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return dnrm2(&n, &x[0], &incx)

cpdef double _test_dzasum(double complex[:] zx) noexcept nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return dzasum(&n, &zx[0], &incx)

cpdef double _test_dznrm2(double complex[:] x) noexcept nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return dznrm2(&n, &x[0], &incx)

cpdef int _test_icamax(float complex[:] cx) noexcept nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return icamax(&n, &cx[0], &incx)

cpdef int _test_idamax(double[:] dx) noexcept nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return idamax(&n, &dx[0], &incx)

cpdef int _test_isamax(float[:] sx) noexcept nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
    return isamax(&n, &sx[0], &incx)

cpdef int _test_izamax(double complex[:] zx) noexcept nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return izamax(&n, &zx[0], &incx)

cpdef float _test_sasum(float[:] sx) noexcept nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
    return sasum(&n, &sx[0], &incx)

cpdef float _test_scasum(float complex[:] cx) noexcept nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return scasum(&n, &cx[0], &incx)

cpdef float _test_scnrm2(float complex[:] x) noexcept nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return scnrm2(&n, &x[0], &incx)

cpdef float _test_sdot(float[:] sx, float[:] sy) noexcept nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
        int incy = sy.strides[0] // sizeof(sy[0])
    return sdot(&n, &sx[0], &incx, &sy[0], &incy)

cpdef float _test_snrm2(float[:] x) noexcept nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return snrm2(&n, &x[0], &incx)

cpdef double complex _test_zdotc(double complex[:] zx, double complex[:] zy) noexcept nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return zdotc(&n, &zx[0], &incx, &zy[0], &incy)

cpdef double complex _test_zdotu(double complex[:] zx, double complex[:] zy) noexcept nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return zdotu(&n, &zx[0], &incx, &zy[0], &incy)
"""

lapack_py_wrappers = """

# Python accessible wrappers for testing:

def _test_dlamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return dlamch(cmach_char)

def _test_slamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return slamch(cmach_char)

cpdef double complex _test_zladiv(double complex zx, double complex zy) noexcept nogil:
    return zladiv(&zx, &zy)

cpdef float complex _test_cladiv(float complex cx, float complex cy) noexcept nogil:
    return cladiv(&cx, &cy)
"""


def generate_pyx(sigs, header_name, accelerate):
    """Generate pyx file with BLAS/LAPACK func declarations and tests."""
    decls = "\n".join([
        pyx_decl(**sig, header_name=header_name, accelerate=accelerate)
        for sig in sigs])
    names = "\n- ".join([sig['name'] for sig in sigs])
    if 'blas' in header_name:
        return blas_pyx_preamble.format(names) + decls + blas_py_wrappers
    return lapack_pyx_preamble.format(names) + decls + lapack_py_wrappers


blas_pxd_preamble = """# Within scipy, these wrappers can be used via relative or absolute cimport.
# Examples:
# from ..linalg cimport cython_blas
# from scipy.linalg cimport cython_blas
# cimport scipy.linalg.cython_blas as cython_blas
# cimport ..linalg.cython_blas as cython_blas

# Within SciPy, if BLAS functions are needed in C/C++/Fortran,
# these wrappers should not be used.
# The original libraries should be linked directly.

ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

"""

lapack_pxd_preamble = """# Within SciPy, these wrappers can be used via relative or absolute cimport.
# Examples:
# from ..linalg cimport cython_lapack
# from scipy.linalg cimport cython_lapack
# cimport scipy.linalg.cython_lapack as cython_lapack
# cimport ..linalg.cython_lapack as cython_lapack

# Within SciPy, if LAPACK functions are needed in C/C++/Fortran,
# these wrappers should not be used.
# The original libraries should be linked directly.

ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

# Function pointer type declarations for
# gees and gges families of functions.
ctypedef bint cselect1(c*)
ctypedef bint cselect2(c*, c*)
ctypedef bint dselect2(d*, d*)
ctypedef bint dselect3(d*, d*, d*)
ctypedef bint sselect2(s*, s*)
ctypedef bint sselect3(s*, s*, s*)
ctypedef bint zselect1(z*)
ctypedef bint zselect2(z*, z*)

"""


def generate_pxd(sigs, lib_name):
    """Create pxd header file for generated pyx."""
    if lib_name == 'BLAS':
        pxd = [blas_pxd_preamble]
    else:
        pxd = [lapack_pxd_preamble]
    for sig in sigs:
        args = ', '.join([
            ' *'.join(arg) for arg in zip(sig['argtypes'], sig['argnames'])])
        pxd.append(
            f"cdef {sig['return_type']} {sig['name']}({args}) noexcept nogil")
    return '\n'.join(pxd)


def c_header_decl(name, return_type, argnames, argtypes, accelerate):
    """Create C header declarations for Cython to import."""
    fort_macro = 'BLAS_FUNC'
    fort_name = name
    return_type = c_types[return_type]
    args = ', '.join(f'{c_types[t]} *{n}' for t, n in zip(argtypes, argnames))
    if name in wrapped_funcs:
        fort_macro = ''
        name = f'{name}wrp_'
        args = f'{return_type} *ret, ' + args
        return_type = 'void'
    # Use old Accelerate symbols if missing from new API
    elif accelerate:
        if name in ['dcabs1', 'lsame']:
            fort_macro = ''
            fort_name += '_'
        elif name == 'xerbla_array':
            fort_macro = ''
            fort_name += '__'
    return f"{return_type} {fort_macro}({name})({args});\n"


c_preamble = """#ifndef SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#define SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#include "fortran_defs.h"
#include "npy_cblas.h"
"""

lapack_decls = """
typedef int (*_cselect1)(npy_complex64*);
typedef int (*_cselect2)(npy_complex64*, npy_complex64*);
typedef int (*_dselect2)(double*, double*);
typedef int (*_dselect3)(double*, double*, double*);
typedef int (*_sselect2)(float*, float*);
typedef int (*_sselect3)(float*, float*, float*);
typedef int (*_zselect1)(npy_complex128*);
typedef int (*_zselect2)(npy_complex128*, npy_complex128*);
"""

cpp_guard_begin = """
#ifdef __cplusplus
extern "C" {
#endif

"""

cpp_guard_end = """
#ifdef __cplusplus
}
#endif
#endif
"""


def generate_c_header(sigs, lib_name, accelerate):
    """Generate complete C header file with declarations and includes."""
    header = [c_preamble.format(lib=lib_name)]
    if lib_name == 'LAPACK':
        header.append(lapack_decls)
    header.append(cpp_guard_begin)
    header += [c_header_decl(**sig, accelerate=accelerate) for sig in sigs]
    header.append(cpp_guard_end)
    return "".join(header)


def newer(source, target):
    """
    Return true if 'source' exists and is more recently modified than
    'target', or if 'source' exists and 'target' doesn't.  Return false if
    both exist and 'target' is the same age or younger than 'source'.
    """
    if not os.path.exists(source):
        raise ValueError("file '%s' does not exist" % os.path.abspath(source))
    if not os.path.exists(target):
        return 1

    mtime1 = os.stat(source)[ST_MTIME]
    mtime2 = os.stat(target)[ST_MTIME]

    return mtime1 > mtime2


def all_newer(src_files, dst_files):
    return all(os.path.exists(dst) and newer(dst, src)
               for dst in dst_files for src in src_files)


def make_all(outdir,
             blas_signature_file="cython_blas_signatures.txt",
             lapack_signature_file="cython_lapack_signatures.txt",
             blas_name="cython_blas",
             lapack_name="cython_lapack",
             blas_header_name="_blas_subroutines.h",
             lapack_header_name="_lapack_subroutines.h",
             accelerate=False):

    src_files = (os.path.abspath(__file__),
                 blas_signature_file,
                 lapack_signature_file)
    dst_files = (blas_name + '.pyx',
                 blas_name + '.pxd',
                 blas_header_name,
                 lapack_name + '.pyx',
                 lapack_name + '.pxd',
                 lapack_header_name)
    dst_files = (os.path.join(outdir, f) for f in dst_files)

    os.chdir(BASE_DIR)

    if all_newer(src_files, dst_files):
        print("scipy/linalg/_generate_pyx.py: all files up-to-date")
        return

    comments = ["This file was generated by _generate_pyx.py.\n",
                "Do not edit this file directly.\n"]
    ccomment = ''.join(['/* ' + line.rstrip() + ' */\n'
                        for line in comments]) + '\n'
    pyxcomment = ''.join(['# ' + line for line in comments]) + '\n'

    with open(blas_signature_file) as f:
        blas_sigs = f.readlines()
    blas_sigs = read_signatures(blas_sigs)
    blas_pyx = generate_pyx(blas_sigs, blas_header_name, accelerate)
    with open(os.path.join(outdir, blas_name + '.pyx'), 'w') as f:
        f.write(pyxcomment)
        f.write(blas_pyx)
    blas_pxd = generate_pxd(blas_sigs, 'BLAS')
    with open(os.path.join(outdir, blas_name + '.pxd'), 'w') as f:
        f.write(pyxcomment)
        f.write(blas_pxd)
    blas_c_header = generate_c_header(blas_sigs, 'BLAS', accelerate)
    with open(os.path.join(outdir, blas_header_name), 'w') as f:
        f.write(ccomment)
        f.write(blas_c_header)

    with open(lapack_signature_file) as f:
        lapack_sigs = f.readlines()
    lapack_sigs = read_signatures(lapack_sigs)
    lapack_pyx = generate_pyx(lapack_sigs, lapack_header_name, accelerate)
    with open(os.path.join(outdir, lapack_name + '.pyx'), 'w') as f:
        f.write(pyxcomment)
        f.write(lapack_pyx)
    lapack_pxd = generate_pxd(lapack_sigs, 'LAPACK')
    with open(os.path.join(outdir, lapack_name + '.pxd'), 'w') as f:
        f.write(pyxcomment)
        f.write(lapack_pxd)
    lapack_c_header = generate_c_header(lapack_sigs, 'LAPACK', accelerate)
    with open(os.path.join(outdir, lapack_header_name), 'w') as f:
        f.write(ccomment)
        f.write(lapack_c_header)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    parser.add_argument("-a", "--accelerate", action="store_true",
                        help="Whether to use new Accelerate (macOS 13.3+)")
    args = parser.parse_args()

    if not args.outdir:
        raise ValueError("Missing `--outdir` argument to _generate_pyx.py")
    else:
        outdir_abs = os.path.join(os.getcwd(), args.outdir)

    make_all(outdir_abs, accelerate=args.accelerate)
