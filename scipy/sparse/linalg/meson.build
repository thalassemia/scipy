py3.install_sources([
    '__init__.py',
    '_expm_multiply.py',
    '_interface.py',
    '_matfuncs.py',
    '_norm.py',
    '_onenormest.py',
    '_svdp.py',
    'dsolve.py',
    'eigen.py',
    'interface.py',
    'isolve.py',
    'matfuncs.py'
  ],
  subdir: 'scipy/sparse/linalg'
)

blas_wrapper = custom_target('blas_wrapper',
  output: [
    '_blas_subroutines.c',
    '_lapack_subroutines.c'
  ],
  input: '_generate_blas_wrapper.py',
  command: [py3, '@INPUT@', '-o', '@OUTDIR@', '-s', blas_symbol_suffix],
  depend_files: [
    '../../linalg/cython_blas_signatures.txt',
    '../../linalg/cython_lapack_signatures.txt',
  ]
)

cwrappers = static_library('cwrappers',
  [blas_wrapper],
  dependencies: [blas, lapack, np_dep, py3_dep],
)

subdir('_propack')
subdir('_isolve')
subdir('_dsolve')
subdir('_eigen')
subdir('tests')
