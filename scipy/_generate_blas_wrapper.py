import argparse
from operator import itemgetter
import os

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
LINALG_DIR = os.path.abspath(os.path.join(CURR_DIR, "linalg"))

c_types = {'int': 'F_INT',
           'c': 'npy_complex64',
           'd': 'double',
           's': 'float',
           'z': 'npy_complex128',
           'char': 'char',
           'bint': 'int',
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


def filter_lines(lines, funcs):
    lines = [line for line in map(str.strip, lines)
                      if line and not line.startswith('#')]
    func_sigs = [split_signature(line) for line in lines
                                           if line.split(' ')[0] != 'void']
    if funcs:
        func_sigs = [f for f in func_sigs if f[0] in funcs]
    sub_sigs = [split_signature(line) for line in lines
                                          if line.split(' ')[0] == 'void']
    if funcs:
        sub_sigs = [f for f in sub_sigs if f[0] in funcs]
    all_sigs = list(sorted(func_sigs + sub_sigs, key=itemgetter(0)))
    return func_sigs, sub_sigs, all_sigs


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


wrapped_c_func_template = """
void {fort_macro}({fort_name})({return_type} *ret, {args});
{return_type} F_FUNC({name}, {upname})({args}){{
    {return_type} ret;
    {fort_macro}({fort_name})(&ret,{f_args});
    return ret;
}}
"""


def c_func_decl(name, return_type, args, suffix):
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
            def_args = ','.join(['int *n', args])
            extra_setup = 'int n = 1;'
            f_args = ','.join(['&n', f_args])
        elif name in ['cdotc', 'cdotu', 'zdotc', 'zdotu', 'cladiv', 'zladiv']:
            name = 'w' + name
            return wrapped_c_func_template.format(name=name, upname=name.upper(),
                                                return_type=return_type, args=args,
                                                f_args=f_args, fort_name=fort_name,
                                                fort_macro=fort_macro)
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
#ifdef HAVE_BLAS_ILP64
#define F_INT npy_int64
#else
#define F_INT int
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


def generate_c_header(func_sigs, sub_sigs, all_sigs, lib_name, suffix):
    funcs = "".join(c_func_decl(sig[0], sig[1], sig[2], suffix) for sig in func_sigs)
    subs = "\n" + "".join(c_sub_decl(sig[0], sig[1], sig[2], suffix) for sig in sub_sigs)
    if lib_name == 'LAPACK':
        preamble = (c_preamble.format(lib=lib_name) + lapack_decls)
    else:
        preamble = c_preamble.format(lib=lib_name)
    return "".join([preamble, cpp_guard, funcs, subs, c_end])


def make_all(outdir,
             blas_signature_file=os.path.join(LINALG_DIR, "cython_blas_signatures.txt"),
             lapack_signature_file=os.path.join(LINALG_DIR, "cython_lapack_signatures.txt"),
             suffix='',
             name='subroutines',
             funcs=None):
    blas_header_name = f'_blas_{name}.c'
    lapack_header_name = f'_lapack_{name}.c'
    comments = ["This file was generated by _generate_pyx.py.\n",
                "Do not edit this file directly.\n"]
    ccomment = ''.join(['/* ' + line.rstrip() + ' */\n'
                        for line in comments]) + '\n'
    with open(blas_signature_file) as f:
        blas_sigs = f.readlines()
    blas_sigs = filter_lines(blas_sigs, funcs)
    blas_c_header = generate_c_header(*(blas_sigs + ('BLAS', suffix)))
    with open(os.path.join(outdir, blas_header_name), 'w') as f:
        f.write(ccomment)
        f.write(blas_c_header)
    with open(lapack_signature_file) as f:
        lapack_sigs = f.readlines()
    lapack_sigs = filter_lines(lapack_sigs, funcs)
    lapack_c_header = generate_c_header(*(lapack_sigs + ('LAPACK', suffix)))
    with open(os.path.join(outdir, lapack_header_name), 'w') as f:
        f.write(ccomment)
        f.write(lapack_c_header)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    parser.add_argument("-s", "--suffix", type=str,
                        help="Suffix for Apple Accelerate")
    parser.add_argument("-n", "--name", type=str, default='subroutines',
                        help="Suffix for output files")
    parser.add_argument("-f", "--funcs", type=str, nargs="+",
                        help="Names of all funcs to wrap")
    args = parser.parse_args()

    if not args.outdir:
        #raise ValueError(f"Missing `--outdir` argument to _generate_pyx.py")
        # We're dealing with a distutils build here, write in-place:
        outdir_abs = os.path.abspath(os.path.dirname(__file__))
    else:
        outdir_abs = os.path.join(os.getcwd(), args.outdir)

    make_all(outdir_abs, suffix=args.suffix, funcs=args.funcs, name=args.name)
