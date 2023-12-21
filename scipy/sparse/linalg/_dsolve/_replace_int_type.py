import argparse
import os
import re

int_re = re.compile(r"\bint\b")
uint_re = re.compile(r"\bunsigned\b \bint\b")

type_include = """
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER) && (_MSC_VER < 1900)
#include "msc_stdint.h"
#else
#include <stdint.h>
#endif
/* Define my integer type int_t */
#ifdef HAVE_BLAS_ILP64
#define NPY_INT_TYPE NPY_INT64
typedef int64_t int_t;
typedef uint64_t uint_t;
#else
#define NPY_INT_TYPE NPY_INT
typedef int int_t;
typedef unsigned int uint_t;
#endif
"""

def replace_types(src_files, out_dir):
    for src_file in src_files:
        code = []
        with open(src_file, "r") as f:
            code += f.readlines()
        code = "".join(code)
        code = uint_re.sub(lambda mobj: 'uint_t', code)
        code = int_re.sub(lambda mobj: 'int_t', code)
        code = "".join([type_include, code])
        with open(os.path.join(out_dir, os.path.split(src_file)[1]), "w") as f:
            f.writelines(code)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", type=str, nargs="+",
                        help="Files to do string substitutions on.")
    parser.add_argument("-o", "--out-dir", type=str,
                        help="Directory to write output files to.")
    args = parser.parse_args()
    replace_types(args.files, args.out_dir)
