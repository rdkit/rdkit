#!/bin/bash

# create binary *.deb packages

#set -x # DEBUG

mkdir -p build
cd build
cmake \
    -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_BUILD_AVALON_SUPPORT=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DRDK_BUILD_CAIRO_SUPPORT=ON \
    -DCMAKE_INSTALL_PREFIX=/usr \
    ../
nprocs=`getconf _NPROCESSORS_ONLN`
make -j $nprocs
cpack -G DEB

# # to install all necessary dependencies on Ubuntu
#sudo apt install -y \
#  curl \
#  wget \
#  libboost-all-dev \
#  libcairo2-dev \
#  cmake \
#  git \
#  g++ \
#  libeigen3-dev \
#  python3 \
#  libpython3-all-dev \
#  python3-numpy \
#  python3-pip \
#  python3-pil \
#  python3-six \
#  python3-pandas

# # to install the freshly built rdkit packages
# sudo dpkg -i *.deb
