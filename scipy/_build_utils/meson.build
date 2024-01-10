if generate_blas_wrappers
  blas_sigs = fs.read('../linalg/cython_blas_signatures.txt').split('\n')
  lapack_sigs = fs.read('../linalg/cython_lapack_signatures.txt').split('\n')
  blas_lapack_names = []
  foreach sig : blas_sigs + lapack_sigs
    if sig.startswith('#') or sig == ''
      continue
    endif
    name = sig.split('(')[0].split(' ')[1]
    blas_lapack_names += name + '.c'
  endforeach

  blas_wrapper_gen = files('_generate_blas_wrapper.py')
  blas_lapack_wrappers = custom_target('blas_lapack_wrappers',
      output: [blas_lapack_names],
      command: [py3, blas_wrapper_gen, '-o', '@OUTDIR@', '-s', blas_symbol_suffix],
      depend_files: [
        '../linalg/cython_blas_signatures.txt',
        '../linalg/cython_lapack_signatures.txt',
      ]
  )

  blas_lapack_wrapper_lib = static_library('blas_lapack_wrapper_lib',
    blas_lapack_wrappers,
    c_args: numpy_nodepr_api,
    dependencies: [blas_dep, lapack_dep, np_dep, py3_dep],
    include_directories: ['../linalg/', './src'],
  )
else
  blas_lapack_wrapper_lib = []
endif

blas_complex_wrp = declare_dependency(
  sources: ['src/complex_blas.c'],
  dependencies: [blas_dep, lapack_dep, np_dep, py3_dep],
  include_directories: ['../linalg/', './src'],
)
