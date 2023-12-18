arpack_sources = [
  'ARPACK/SRC/cgetv0.f',
  'ARPACK/SRC/cnaitr.f',
  'ARPACK/SRC/cnapps.f',
  'ARPACK/SRC/cnaup2.f',
  'ARPACK/SRC/cnaupd.f',
  'ARPACK/SRC/cneigh.f',
  'ARPACK/SRC/cneupd.f',
  'ARPACK/SRC/cngets.f',
  'ARPACK/SRC/csortc.f',
  'ARPACK/SRC/cstatn.f',
  'ARPACK/SRC/dgetv0.f',
  'ARPACK/SRC/dnaitr.f',
  'ARPACK/SRC/dnapps.f',
  'ARPACK/SRC/dnaup2.f',
  'ARPACK/SRC/dnaupd.f',
  'ARPACK/SRC/dnconv.f',
  'ARPACK/SRC/dneigh.f',
  'ARPACK/SRC/dneupd.f',
  'ARPACK/SRC/dngets.f',
  'ARPACK/SRC/dsaitr.f',
  'ARPACK/SRC/dsapps.f',
  'ARPACK/SRC/dsaup2.f',
  'ARPACK/SRC/dsaupd.f',
  'ARPACK/SRC/dsconv.f',
  'ARPACK/SRC/dseigt.f',
  'ARPACK/SRC/dsesrt.f',
  'ARPACK/SRC/dseupd.f',
  'ARPACK/SRC/dsgets.f',
  'ARPACK/SRC/dsortc.f',
  'ARPACK/SRC/dsortr.f',
  'ARPACK/SRC/dstatn.f',
  'ARPACK/SRC/dstats.f',
  'ARPACK/SRC/dstqrb.f',
  'ARPACK/SRC/sgetv0.f',
  'ARPACK/SRC/snaitr.f',
  'ARPACK/SRC/snapps.f',
  'ARPACK/SRC/snaup2.f',
  'ARPACK/SRC/snaupd.f',
  'ARPACK/SRC/snconv.f',
  'ARPACK/SRC/sneigh.f',
  'ARPACK/SRC/sneupd.f',
  'ARPACK/SRC/sngets.f',
  'ARPACK/SRC/ssaitr.f',
  'ARPACK/SRC/ssapps.f',
  'ARPACK/SRC/ssaup2.f',
  'ARPACK/SRC/ssaupd.f',
  'ARPACK/SRC/ssconv.f',
  'ARPACK/SRC/sseigt.f',
  'ARPACK/SRC/ssesrt.f',
  'ARPACK/SRC/sseupd.f',
  'ARPACK/SRC/ssgets.f',
  'ARPACK/SRC/ssortc.f',
  'ARPACK/SRC/ssortr.f',
  'ARPACK/SRC/sstatn.f',
  'ARPACK/SRC/sstats.f',
  'ARPACK/SRC/sstqrb.f',
  'ARPACK/SRC/zgetv0.f',
  'ARPACK/SRC/znaitr.f',
  'ARPACK/SRC/znapps.f',
  'ARPACK/SRC/znaup2.f',
  'ARPACK/SRC/znaupd.f',
  'ARPACK/SRC/zneigh.f',
  'ARPACK/SRC/zneupd.f',
  'ARPACK/SRC/zngets.f',
  'ARPACK/SRC/zsortc.f',
  'ARPACK/SRC/zstatn.f',
  'ARPACK/UTIL/cmout.f',
  'ARPACK/UTIL/cvout.f',
  'ARPACK/UTIL/dmout.f',
  'ARPACK/UTIL/dvout.f',
  'ARPACK/UTIL/icnteq.f',
  'ARPACK/UTIL/icopy.f',
  'ARPACK/UTIL/iset.f',
  'ARPACK/UTIL/iswap.f',
  'ARPACK/UTIL/ivout.f',
  'ARPACK/UTIL/second_NONE.f',
  'ARPACK/UTIL/smout.f',
  'ARPACK/UTIL/svout.f',
  'ARPACK/UTIL/zmout.f',
  'ARPACK/UTIL/zvout.f'
]

# Building ARPACK yields a ton of rank mismatch (scalar and rank-1) warnings
# that cannot be suppressed with a more specific flag.
_suppress_all_warnings = ff.get_supported_arguments('-w')

arpack_lib = static_library('arpack_lib',
  arpack_sources,
  fortran_args: [fortran_ignore_warnings, _suppress_all_warnings, _fortran_64_bit_int],
  dependencies: [lapack, blas],
  include_directories: ['ARPACK/SRC']
)

arpack_blas_wrapper = custom_target('arpack_blas_wrapper',
  output: [
    '_blas_arpack.c',
    '_lapack_arpack.c'
  ],
  command: [py3, generate_blas_wrapper, '-f', 'caxpy', 'ccopy', 'cgemv', 'cgeqr2', 'cgeru', 'clacpy', 'clahqr', 'clanhs', 'clarnv', 'clartg', 'clascl', 'claset', 'cscal', 'csscal', 'ctrevc', 'ctrmm', 'ctrsen', 'cunm2r', 'daxpy', 'dcopy', 'ddot', 'dgemv', 'dgeqr2', 'dger', 'dlabad', 'dlacpy', 'dlaev2', 'dlahqr', 'dlamch', 'dlanhs', 'dlanst', 'dlapy2', 'dlarf', 'dlarfg', 'dlarnv', 'dlartg', 'dlascl', 'dlaset', 'dlasr', 'dnrm2', 'dorm2r', 'dscal', 'dsteqr', 'dswap', 'dtrevc', 'dtrmm', 'dtrsen', 'dznrm2', 'saxpy', 'scnrm2', 'scopy', 'sdot', 'sgemv', 'sgeqr2', 'sger', 'slabad', 'slacpy', 'slaev2', 'slahqr', 'slamch', 'slanhs', 'slanst', 'slapy2', 'slarf', 'slarfg', 'slarnv', 'slartg', 'slascl', 'slaset', 'slasr', 'snrm2', 'sorm2r', 'sscal', 'ssteqr', 'sswap', 'strevc', 'strmm', 'strsen', 'cdotc', 'zdotc', 'zaxpy', 'zcopy', 'zdscal', 'zgemv', 'zgeqr2', 'zgeru', 'zlacpy', 'zlahqr', 'zlanhs', 'zlarnv', 'zlartg', 'zlascl', 'zlaset', 'zscal', 'ztrevc', 'ztrmm', 'ztrsen', 'zunm2r', '-o', '@OUTDIR@', '-s', blas_symbol_suffix, '-n', 'arpack'],
  depend_files: [
    '../../../../linalg/cython_blas_signatures.txt',
    '../../../../linalg/cython_lapack_signatures.txt',
  ]
)

arpack_module = custom_target('arpack_module',
  output: ['_arpackmodule.c', '_arpack-f2pywrappers.f'],
  input: 'arpack.pyf.src',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_arpack = py3.extension_module('_arpack',
  [arpack_module, arpack_blas_wrapper],
  c_args: numpy_nodepr_api,
  link_with: [arpack_lib],
  link_args: version_link_args,
  dependencies: [lapack, blas, fortranobject_dep, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/sparse/linalg/_eigen/arpack',
  include_directories: ['../../../../linalg']
)

py3.install_sources([
    '__init__.py',
    'arpack.py',
    'ARPACK/COPYING'
  ],
subdir: 'scipy/sparse/linalg/_eigen/arpack'
)

subdir('tests')
