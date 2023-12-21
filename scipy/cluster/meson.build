slib = py3.extension_module('_hierarchy',
  linalg_cython_gen.process('_hierarchy.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)

slib = py3.extension_module('_optimal_leaf_ordering',
  linalg_cython_gen.process('_optimal_leaf_ordering.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)

_vq_pyx = custom_target('_vq',
  output: '_vq.pyx',
  input: '_vq.pyx.in',
  command: [py3, tempita, '-o', '@OUTDIR@', '@INPUT@', 'i_type=' + blas_int_type]
)

slib = py3.extension_module('_vq',
  linalg_cython_gen.process(_vq_pyx),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster'
)


py3.install_sources([
    '__init__.py',
    'hierarchy.py',
    'vq.py'
  ],
  subdir: 'scipy/cluster'
)

subdir('tests')
