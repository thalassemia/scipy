set -xe


printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
echo PKG_CONFIG_PATH $PKG_CONFIG_PATH
PKG_CONFIG_PATH=$PROJECT_DIR/.openblas
rm -rf $PKG_CONFIG_PATH
mkdir -p $PKG_CONFIG_PATH
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PKG_CONFIG_PATH/scipy-openblas.pc
# Copy the shared objects to a path under $PKG_CONFIG_PATH, the build
# will point $LD_LIBRARY_PATH there and then auditwheel/delocate-wheel will
# pull these into the wheel. Use python to avoid windows/posix problems
python <<EOF
import os, scipy_openblas32, shutil
srcdir = os.path.join(os.path.dirname(scipy_openblas32.__file__), "lib")
shutil.copytree(srcdir, os.path.join("$PKG_CONFIG_PATH", "lib"))
srcdir = os.path.join(os.path.dirname(scipy_openblas32.__file__), ".dylibs")
if os.path.exists(srcdir):  # macosx delocate
    shutil.copytree(srcdir, os.path.join("$PKG_CONFIG_PATH", ".dylibs"))
EOF


