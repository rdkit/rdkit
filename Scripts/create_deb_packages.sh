#!/bin/bash

# create binary *.deb packages

#set -x # DEBUG

mkdir -p build
cd build
cmake -Wno-dev \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_BUILD_AVALON_SUPPORT=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DCMAKE_INSTALL_PREFIX=/usr \
    ../
nprocs=`getconf _NPROCESSORS_ONLN`
make -j $nprocs
cpack -G DEB

# # to install all necessary dependencies on Ubuntu
# sudo apt-get install \
#      fonts-freefont-ttf \
#      libboost-python1.58.0 \
#      libboost-regex1.58.0 \
#      libboost-system1.58.0 \
#      libboost-thread1.58.0 \
#      libc6 \
#      libgcc1 \
#      libpython2.7 \
#      libstdc++6 \
#      python
# # to install the freshly built rdkit packages
# sudo dpkg -i *.deb
