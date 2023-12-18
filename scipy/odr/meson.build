odrpack = static_library('odrpack',
  [
    'odrpack/d_lpk.f',
    'odrpack/d_mprec.f',
    'odrpack/d_odr.f',
    'odrpack/dlunoc.f'
  ],
  fortran_args: [_fflag_Wno_conversion, _fortran_64_bit_int],  # silence "conversion from REAL(8) to INTEGER(4)"
)

odrpack_blas_wrapper = custom_target('odrpack_blas_wrapper',
  output: [
    '_blas_odrpack.c',
    '_lapack_odrpack.c'
  ],
  command: [py3, generate_blas_wrapper, '-f', 'dasum', 'daxpy', 'dcopy', 'ddot', 'dnrm2', 'drot', 'drotg', 'dscal', 'dswap', 'idamax', '-o', '@OUTDIR@', '-s', blas_symbol_suffix, '-n', 'odrpack'],
  depend_files: [
    '../linalg/cython_blas_signatures.txt',
    '../linalg/cython_lapack_signatures.txt',
  ]
)

py3.extension_module('__odrpack',
  ['__odrpack.c', odrpack_blas_wrapper],
  link_with: odrpack,
  link_args: version_link_args,
  c_args: numpy_nodepr_api,
  dependencies: [blas, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/odr',
  include_directories: ['../linalg']
)


python_sources = [
  '__init__.py',
  '_add_newdocs.py',
  '_models.py',
  '_odrpack.py',
  'models.py',
  'odrpack.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/odr'
)

subdir('tests')
