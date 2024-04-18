set -xe

PROJECT_DIR="$1"

# python $PROJECT_DIR/tools/wheels/check_license.py
if [[ $(uname) == "Linux" ]] ; then
    python $PROJECT_DIR/tools/openblas_support.py --check_version
# Set in .github/workflows/wheels.yml for x86 Accelerate wheel
# elif [[ $VTOOL_PATCH_MIN_VER ]]
#     SCIPY_LIB=$(python -c "import os; print(os.path.dirname(os.__file__) + '/scipy')")
#     # Required to run test suite on macOS 13 runner
#     find $SCIPY_LIB -type f -name '*.so' -exec vtool -set-version-min macos 13.6 14.2 -o "{}" "{}" \;
fi
echo $?

python -c "import sys; import scipy; sys.exit(not scipy.test())"
echo $?
