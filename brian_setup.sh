set -e
conda deactivate || /usr/bin/true
echo "variables"
export compiler="clangxx_osx-arm64=18.1"
export compiler_version=16
export boost_version=1.84
export number_of_cores=`sysctl -n hw.ncpu`
export target_platform=15.1
export python="python=3.11"

echo "conda create"
export CONDA=/Users/guillaume-osmo/miniconda3

if [ ! -d $CONDA ] ; then
   #echo "making conda"
   # wget https://github.com/conda-forge/miniforge/releases/download/24.9.0-0/Miniforge3-24.9.0-0-MacOSX-arm64.sh
   # bash Miniforge3-24.9.0-0-MacOSX-arm64.sh -b -p ${CONDA}

   #source ${CONDA}/etc/profile.d/conda.sh
   #conda update -q conda
   conda config --set solver libmamba
   conda config --set channel_priority strict
   conda config --add channels conda-forge

   echo "conda create -y --name rdkit_build_brian $python $compiler libcxx cmake \
        libboost=$boost_version libboost-devel=$boost_version \
        libboost-python=$boost_version libboost-python-devel=$boost_version \
        qt \
        numpy matplotlib=3.8 cairo pillow eigen pandas=2.1  \
        jupyter=1.0 ipython=8.20 sphinx myst-parser pytest nbval"

   conda create --yes --name rdkit_build_brian $python $compiler libcxx cmake \
         libboost=$boost_version libboost-devel=$boost_version \
         libboost-python=$boost_version libboost-python-devel=$boost_version \
         qt \
         numpy matplotlib=3.8 cairo pillow eigen pandas=2.1  \
         jupyter=1.0 ipython=8.20 sphinx myst-parser pytest nbval

fi

echo "activate"
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX$target_platform.sdk
export CONDA_BUILD_SYSROOT=${SDKROOT}
export CXXFLAGS="${CXXFLAGS} -Wall -Wextra -Werror"

conda activate rdkit_build
conda install -c conda-forge openblas
echo "mac build"
rm -rf mac-build
mkdir mac-build && cd mac-build && \
    cmake .. \
          -DCMAKE_BUILD_TYPE=Debug \
          -DRDK_INSTALL_INTREE=ON \
          -DRDK_INSTALL_STATIC_LIBS=OFF \
          -DRDK_BUILD_CPP_TESTS=ON \
          -DRDK_BUILD_PYTHON_WRAPPERS=ON \
          -DRDK_BUILD_COORDGEN_SUPPORT=ON \
          -DRDK_BUILD_MAEPARSER_SUPPORT=ON \
          -DRDK_OPTIMIZE_POPCNT=ON \
          -DRDK_BUILD_TEST_GZIP=ON \
          -DRDK_BUILD_FREESASA_SUPPORT=ON \
          -DRDK_BUILD_AVALON_SUPPORT=ON \
          -DRDK_BUILD_INCHI_SUPPORT=ON \
          -DRDK_BUILD_YAEHMOP_SUPPORT=ON \
          -DRDK_BUILD_XYZ2MOL_SUPPORT=ON \
          -DRDK_BUILD_CAIRO_SUPPORT=ON \
          -DRDK_BUILD_QT_SUPPORT=ON \
          -DRDK_BUILD_SWIG_WRAPPERS=OFF \
          -DRDK_SWIG_STATIC=OFF \
          -DRDK_BUILD_THREADSAFE_SSS=ON \
          -DRDK_TEST_MULTITHREADED=ON \
          -DRDK_BUILD_CFFI_LIB=ON \
          -DRDK_BUILD_OSMORDRED_SUPPORT=ON \
          -DCMAKE_OSX_SYSROOT=${SDKROOT} \
          -DCMAKE_OSX_DEPLOYMENT_TARGET=$target_platform