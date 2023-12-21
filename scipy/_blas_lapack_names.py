import argparse
from operator import itemgetter
import os

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
LINALG_DIR = os.path.abspath(os.path.join(CURR_DIR, "linalg"))

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
    all_sigs = list(sorted(func_sigs + sub_sigs, key=itemgetter(0)))
    return func_sigs, sub_sigs, all_sigs

def make_all(outdir,
             blas_signature_file=os.path.join(LINALG_DIR, "cython_blas_signatures.txt"),
             lapack_signature_file=os.path.join(LINALG_DIR, "cython_lapack_signatures.txt")):
    with open(blas_signature_file) as f:
        blas_sigs = f.readlines()
    _, _, blas_sigs = filter_lines(blas_sigs)
    with open(lapack_signature_file) as f:
        lapack_sigs = f.readlines()
    _, _, lapack_sigs = filter_lines(lapack_sigs)
    func_names = []
    for sig in blas_sigs + lapack_sigs:
        func_names.append(f'_{sig[0]}.c')
    with open(os.path.join(outdir, 'blas_lapack_func_names.txt'), 'w') as f:
        f.writelines("\n".join(func_names))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    args = parser.parse_args()

    if not args.outdir:
        #raise ValueError(f"Missing `--outdir` argument to _generate_pyx.py")
        # We're dealing with a distutils build here, write in-place:
        outdir_abs = os.path.abspath(os.path.dirname(__file__))
    else:
        outdir_abs = os.path.join(os.getcwd(), args.outdir)

    make_all(outdir_abs)
