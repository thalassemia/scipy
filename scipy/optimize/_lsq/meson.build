givens_elimination_pyx = custom_target('givens_elimination',
  output: 'givens_elimination.pyx',
  input: 'givens_elimination.pyx.in',
  command: [py3, tempita, '-o', '@OUTDIR@', '@INPUT@', 'i_type=' + blas_int_type]
)

_lsq = py3.extension_module('givens_elimination',
  linalg_cython_gen.process(givens_elimination_pyx),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/optimize/_lsq'
)

py3.install_sources([
    '__init__.py',
    'bvls.py',
    'common.py',
    'dogbox.py',
    'least_squares.py',
    'lsq_linear.py',
    'trf.py',
    'trf_linear.py'
  ],
  subdir: 'scipy/optimize/_lsq'
)
