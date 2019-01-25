#!/bin/bash

#set -x # DEBUG

# launch the unit tests suite.
# run this script like this: ./Scripts/test.sh

RDBASE=$PWD
PYTHONPATH=$RDBASE
LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH

# run them in parallel
nprocs=`getconf _NPROCESSORS_ONLN`

cd build && ctest --output-on-failure -j ${nprocs}
