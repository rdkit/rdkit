#!/bin/bash

# create *.deb packages for Ubuntu

#set -x # DEBUG

mkdir -p build
cd build
cmake -DRDK_INSTALL_INTREE=OFF -DCMAKE_INSTALL_PREFIX=/usr ../
nprocs=`getconf _NPROCESSORS_ONLN`
make -j $nprocs
cpack -G DEB

# # to install them
# sudo dpkg -i *.deb
