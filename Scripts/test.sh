#!/bin/bash

#set -x # DEBUG

# launch the unit tests suite.
# run this script like this: ./Scripts/test.sh [build_dir]
# if no build_dir is specified, ./build is assumed to be the default
# build directory

build_dir="build"
if [ "$#" -eq 1 ]; then
    build_dir=$1
fi

RDBASE=$PWD
PYTHONPATH=$RDBASE
LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH

# run tests in parallel
nprocs=`getconf _NPROCESSORS_ONLN`

cd $build_dir && ctest --output-on-failure -j ${nprocs}
