import argparse
from operator import itemgetter
import os

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
LINALG_DIR = os.path.abspath(os.path.join(CURR_DIR, "..", "linalg"))

c_types = {'int': 'F_INT',
           'c': 'npy_complex64',
           'd': 'double',
           's': 'float',
           'z': 'npy_complex128',
           'char': 'char',
           'bint': 'F_INT',
           'cselect1': '_cselect1',
           'cselect2': '_cselect2',
           'dselect2': '_dselect2',
           'dselect3': '_dselect3',
           'sselect2': '_sselect2',
           'sselect3': '_sselect3',
           'zselect1': '_zselect1',
           'zselect2': '_zselect2'}


def split_signature(sig):
    name_and_type, args = sig[:-1].split('(')
    ret_type, name = name_and_type.split(' ')
    return name, ret_type, args


def filter_lines(lines):
    lines = [line for line in map(str.strip, lines)
                      if line and not line.startswith('#')]
    func_sigs = [split_signature(line) for line in lines
                                           if line.split(' ')[0] != 'void']
    sub_sigs = [split_signature(line) for line in lines
                                          if line.split(' ')[0] == 'void']
    return func_sigs, sub_sigs


def arg_names_and_types(args):
    return zip(*[arg.split(' *') for arg in args.split(', ')])


def make_c_args(args):
    types, names = arg_names_and_types(args)
    types = [c_types[arg] for arg in types]
    return ', '.join(f'{t} *{n}' for t, n in zip(types, names)), ','.join(names)


c_func_template = """
{return_type} {fort_macro}({fort_name})({def_args});
{return_type} F_FUNC({name}, {upname})({args}){{
    {extra_setup}
    return {fort_macro}({fort_name})({f_args});
}}
"""


complex_dot_template = """
void CBLAS_FUNC({fort_name})({f_args}, {return_type} *ret);
{return_type} F_FUNC({name}, {upname})({args}){{
    {return_type} ret;
    CBLAS_FUNC({fort_name})({c_args},&ret);
    return ret;
}}
"""


ladiv_template = """
void BLAS_FUNC({fort_name})({f_args});
void F_FUNC({name}, {upname})({return_type} *ret, {args}){{
    {float_type} zr, zi, xr, xi, yr, yi;
    xr = npy_creal{f}(*x);
    xi = npy_cimag{f}(*x);
    yr = npy_creal{f}(*y);
    yi = npy_cimag{f}(*y);
    BLAS_FUNC({fort_name})(&xr, &xi, &yr, &yi, &zr, &zi);
    NPY_CSETREAL(ret, zr);
    NPY_CSETIMAG(ret, zi);
}}
"""


def c_func_decl(name, return_type, args, suffix, g77):
    orig_args = args
    args, f_args = make_c_args(args)
    def_args = args
    return_type = c_types[return_type]
    fort_name = name
    fort_macro = 'BLAS_FUNC'
    extra_setup = ''
    if '$NEWLAPACK' in suffix:
        if name == 'dcabs1':
            fort_macro = ''
        elif name == 'lsame':
            fort_name = 'lsamen'
            def_args = ','.join(['F_INT *n', args])
            extra_setup = 'F_INT n = 1;'
            f_args = ','.join(['&n', f_args])
    if g77 and name in ['cdotc', 'cdotu', 'zdotc', 'zdotu']:
        fort_name = f'cblas_{name}_sub'
        name = 'w' + name
        types, names = arg_names_and_types(orig_args)
        types = [c_types[t] for t in types]
        f_args = f'{types[0]} {names[0]}, {types[1]} *{names[1]}, {types[2]} {names[2]}, {types[3]} *{names[3]}, {types[4]} {names[4]}'
        c_args = f'*{names[0]}, {names[1]}, *{names[2]}, {names[3]}, *{names[4]}'
        return complex_dot_template.format(name=name, upname=name.upper(),
            return_type=return_type, args=args, f_args=f_args, 
            fort_name=fort_name, c_args=c_args)
    if g77 and name in ['cladiv', 'zladiv']:
        name = 'w' + name 
        argtypes, argnames = arg_names_and_types(args)
        f_args = []
        if name == 'wcladiv':
            fort_name = 'sladiv'
            float_type = 'float'
            f = 'f'
            f_ret_args = ['float *retr, float * reti']
        else:
            fort_name = 'dladiv'
            float_type = 'double'
            f = ''
            f_ret_args = ['double *retr, double * reti']
        for argtype, argname in zip(argtypes, argnames):
            if argtype == 'npy_complex128':
                f_args.append(f'double *{argname}r, double *{argname}i')
            elif argtype == 'npy_complex64':
                f_args.append(f'float *{argname}r, float *{argname}i')
            else:
                f_args.append(f'{argtype} *{argname}')
        f_args = ', '.join(f_args + f_ret_args)
        return ladiv_template.format(name=name, upname=name.upper(), fort_name=fort_name, 
            f_args=f_args, args=args, float_type=float_type, f=f, 
            return_type=return_type)
    if suffix == '':
        return ''
    return c_func_template.format(name=name, upname=name.upper(),
                                  return_type=return_type, args=args,
                                  f_args=f_args, fort_macro=fort_macro,
                                  fort_name=fort_name, def_args=def_args,
                                  extra_setup=extra_setup)


c_sub_template = """
void {fort_macro}({name})({args});
void F_FUNC({name}, {upname})({args}){{
    {fort_macro}({name})({f_args});
}}
"""


def c_sub_decl(name, return_type, args, suffix):
    # No wrapper required if no suffix
    if suffix == '':
        return ''
    args, f_args = make_c_args(args)
    fort_macro = 'BLAS_FUNC'
    if '$NEWLAPACK' in suffix:
        if name == 'xerbla_array':
            fort_macro = ''
            name += '__'
    return c_sub_template.format(name=name, upname=name.upper(), args=args, fort_macro=fort_macro, f_args=f_args)


c_preamble = """#ifndef SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#define SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#include "fortran_defs.h"
#include "numpy/arrayobject.h"

#include <numpy/npy_math.h>

#ifndef NUMPY_CORE_INCLUDE_NUMPY_NPY_2_COMPLEXCOMPAT_H_
#define NUMPY_CORE_INCLUDE_NUMPY_NPY_2_COMPLEXCOMPAT_H_

#ifndef NPY_CSETREALF
#define NPY_CSETREALF(c, r) (c)->real = (r)
#endif
#ifndef NPY_CSETIMAGF
#define NPY_CSETIMAGF(c, i) (c)->imag = (i)
#endif
#ifndef NPY_CSETREAL
#define NPY_CSETREAL(c, r)  (c)->real = (r)
#endif
#ifndef NPY_CSETIMAG
#define NPY_CSETIMAG(c, i)  (c)->imag = (i)
#endif
#ifndef NPY_CSETREALL
#define NPY_CSETREALL(c, r) (c)->real = (r)
#endif
#ifndef NPY_CSETIMAGL
#define NPY_CSETIMAGL(c, i) (c)->imag = (i)
#endif

#endif
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

cpp_guard = """
#ifdef __cplusplus
extern "C" {
#endif

"""

c_end = """
#ifdef __cplusplus
}
#endif
#endif
"""

comments = ["This file was generated by _generate_pyx.py.\n",
            "Do not edit this file directly.\n"]
ccomment = ''.join(['/* ' + line.rstrip() + ' */\n'
                    for line in comments]) + '\n'


def generate_c_file(func_sigs, sub_sigs, lib_name, suffix, g77, outdir):
    if lib_name == 'LAPACK':
        preamble = (c_preamble.format(lib=lib_name) + lapack_decls)
        out_name = 'lapack_wrappers.c'
    else:
        preamble = c_preamble.format(lib=lib_name)
        out_name = 'blas_wrappers.c'
    funcs_and_subs = [ccomment, preamble, cpp_guard]
    for sig in func_sigs:
        funcs_and_subs.append(c_func_decl(*(sig+(suffix, g77))))
    for sig in sub_sigs:
        funcs_and_subs.append(c_sub_decl(*(sig+(suffix,))))
    funcs_and_subs.append(c_end)
    with open(os.path.join(outdir, out_name), 'w') as func_file:
        func_file.writelines("".join(funcs_and_subs))


def make_all(outdir,
             blas_signature_file=os.path.join(LINALG_DIR, "cython_blas_signatures.txt"),
             lapack_signature_file=os.path.join(LINALG_DIR, "cython_lapack_signatures.txt"),
             suffix='',
             g77=False):
    with open(blas_signature_file) as f:
        blas_sigs = f.readlines()
    blas_sigs = filter_lines(blas_sigs)
    with open(lapack_signature_file) as f:
        lapack_sigs = f.readlines()
    lapack_sigs = filter_lines(lapack_sigs)
    generate_c_file(*(blas_sigs + ('BLAS', suffix, g77, outdir)))
    generate_c_file(*(lapack_sigs + ('LAPACK', suffix, g77, outdir)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    parser.add_argument("-s", "--suffix", type=str,
                        help="Suffix for Apple Accelerate")
    parser.add_argument("-g", "--g77", type=str, default=False,
                        help="Whether to generate g77 wrappers")
    args = parser.parse_args()

    if not args.outdir:
        #raise ValueError(f"Missing `--outdir` argument to _generate_pyx.py")
        # We're dealing with a distutils build here, write in-place:
        outdir_abs = os.path.abspath(os.path.dirname(__file__))
    else:
        outdir_abs = os.path.join(os.getcwd(), args.outdir)

    make_all(outdir_abs, suffix=args.suffix, g77=args.g77)
