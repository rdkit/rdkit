#!/bin/bash

#set -x # DEBUG

# developper build for Linux with Python-3 and conda
# you must have conda in your PATH
# go to https://www.anaconda.com/download/#linux if you don't
# you may need to source ~/.bashrc after the conda install
# so that the default python is the one installed by conda

# configuration
export RDBASE=`pwd`
echo $RDBASE
export PYTHONPATH=$RDBASE

export PYTHON=`which python`
echo $PYTHON
export PY_PREFIX=`$PYTHON -c "import sys; print(sys.prefix)"`
echo $PY_PREFIX
export LD_LIBRARY_PATH=${RDBASE}/lib:${PY_PREFIX}/lib
export PY_SP_DIR=${PY_PREFIX}/lib/python${PYTHON_VERSION}/site-packages
echo $PY_SP_DIR

PYTHON_VERSION=`$PYTHON -c 'import sys; print("%d.%d" % (sys.version_info[0], sys.version_info[1]))'`
# some system-wide dependencies that may be needed
PKGS="gcc make python${PYTHON_VERSION}-dev"
dpkg-query -l $PKGS > /dev/null || sudo apt install $PKGS

conda install -y \
      cmake cairo numpy pillow pandas eigen pkg-config \
      boost-cpp boost py-boost gxx_linux-64

# build
mkdir -p build
cd build
cmake -D RDK_USE_BOOST_REGEX=ON \
      -D PYTHON_EXECUTABLE=$PYTHON \
      -D BOOST_ROOT=$PY_PREFIX \
      -D Boost_NO_SYSTEM_PATHS=ON \
      -D RDK_BUILD_AVALON_SUPPORT=ON \
      -D RDK_BUILD_INCHI_SUPPORT=ON \
      -D RDK_BUILD_THREADSAFE_SSS=on \
      -D RDK_TEST_MULTITHREADED=on \
      -D RDK_INSTALL_STATIC_LIBS=OFF \
      -D RDK_BUILD_SWIG_WRAPPERS=OFF \
      -D RDK_SWIG_STATIC=OFF \
      -D RDK_BUILD_PYTHON_WRAPPERS=ON \
      -D RDK_BUILD_FREESASA_SUPPORT=ON \
      -D PYTHON_NUMPY_INCLUDE_PATH=${PY_SP_DIR}/numpy/core/include \
      ../
nprocs=`getconf _NPROCESSORS_ONLN`
make -j $nprocs

# install
make install

# test
ctest -j $nprocs --output-on-failure
