#!/bin/bash

#set -x # DEBUG

# developper build for Linux with Python-3 and conda

## some system-wide dependencies that may be needed
# sudo apt install gcc make

conda install -y \
      cmake cairo numpy pillow pandas eigen pkg-config \
      boost-cpp boost py-boost gxx_linux-64

PYTHON_VERSION=`python -c 'import sys; print("%d.%d" % (sys.version_info[0], sys.version_info[1]))'`

echo "BEFORE SCRIPT"
export RDBASE=`pwd`
echo $RDBASE
export PYTHONPATH=${RDBASE}
export LD_LIBRARY_PATH=${RDBASE}/lib

export PYTHON=`which python`
echo $PYTHON
export PY_PREFIX=`$PYTHON -c "import sys; print(sys.prefix)"`
echo $PY_PREFIX
export PY_SP_DIR=$PY_PREFIX/lib/python${PYTHON_VERSION}/site-packages
echo $PY_SP_DIR

mkdir -p build
cd build

echo "SCRIPT"
cd $RDBASE
mkdir build
cd build
export REGEX_EXTRA="-DRDK_USE_BOOST_REGEX=ON"
cmake $REGEX_EXTRA \
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
      -D RDK_BUILD_FREESASA_SUPPORT=ON ../
nprocs=`getconf _NPROCESSORS_ONLN`
make -j ${nprocs}
make install
LD_LIBRARY_PATH="$PY_PREFIX/lib:$PREFIX/lib;$SRC_DIR/lib;$LD_LIBRARY_PATH" \
ctest -j ${nprocs} --output-on-failure
