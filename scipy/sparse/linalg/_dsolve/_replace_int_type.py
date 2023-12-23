import argparse
import os
import re

int_re = re.compile(r"\bint\b")
uint_re = re.compile(r"\bunsigned\b \bint\b")
npy_int_re = re.compile(r"\bNPY_INT\b")
int_format_re = re.compile(r"(%.?)d")

type_include = """
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER) && (_MSC_VER < 1900)
#include "msc_stdint.h"
#include "msc_inttypes.h"
#else
#include <stdint.h>
#include <inttypes.h>
#endif
"""

printf_replace = [
    'sutil.c',
    'sreadhb.c',
    'sreadrb.c',
    'sreadtriple.c',
    'smemory.c',
    'sldperm.c',
    'sgstrs.c',
    'cutil.c',
    'creadhb.c',
    'creadrb.c',
    'creadtriple.c',
    'cmemory.c',
    'cldperm.c',
    'cgstrs.c',
    'zutil.c',
    'zreadhb.c',
    'zreadrb.c',
    'zreadtriple.c',
    'zmemory.c',
    'zldperm.c',
    'zgstrs.c',
    'dutil.c',
    'dreadhb.c',
    'dreadrb.c',
    'dreadtriple.c',
    'dmemory.c',
    'dldperm.c',
    'dgstrs.c',
    'util.c',
    'sp_preorder.c',
    'input_error.c',
    'colamd.c'
]

def replace_types(src_files, out_dir, int64):
    for src_file in src_files:
        code = []
        with open(src_file, "r") as f:
            code += f.readlines()
        if int64:
            code = "".join(code)
            code = uint_re.sub(lambda mobj: 'uint64_t', code)
            code = int_re.sub(lambda mobj: 'int64_t', code)
            code = npy_int_re.sub(lambda mobj: 'NPY_INT64', code)
            if os.path.split(src_file)[-1] in printf_replace:
                code = int_format_re.sub(lambda mobj: f'{mobj[1]}" PRId64 "', code)
            code = "".join([type_include, code])
        with open(os.path.join(out_dir, os.path.split(src_file)[1]), "w") as f:
            f.writelines(code)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", type=str, nargs="+",
                        help="Files to do string substitutions on.")
    parser.add_argument("-o", "--out-dir", type=str,
                        help="Directory to write output files to.")
    parser.add_argument("-i", "--int64", action="store_true",
                        help="Use 64-bit integer.")
    args = parser.parse_args()
    replace_types(args.files, args.out_dir, args.int64)
