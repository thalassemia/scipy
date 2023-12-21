generate_blas_wrapper = files('_generate_blas_wrapper.py')
blas_lapack_names = fs.read('../linalg/blas_lapack_func_names.txt').strip().split('\n')
blas_lapack_wrappers = custom_target('blas_lapack_wrappers',
    output: blas_lapack_names,
    command: [py3, generate_blas_wrapper, '-o', '@OUTDIR@', '-s', blas_symbol_suffix, '-g', use_g77_abi.to_string()],
    depend_files: [
      '../linalg/cython_blas_signatures.txt',
      '../linalg/cython_lapack_signatures.txt',
    ]
)
blas_lapack_wrapper_lib = static_library('blas_lapack_wrapper_lib',
  blas_lapack_wrappers,
  c_args: numpy_nodepr_api,
  dependencies: [blas, lapack, np_dep, py3_dep],
  include_directories: ['../linalg/'],
)
