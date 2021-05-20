conda create -n rdkit_dev -c conda-forge -y cmake cairo pillow eigen pkg-config boost-cpp boost py-boost gxx_linux-64 numpy

conda activate rdkit_dev

mkdir build && cd build
cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
    -D PYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy ; print(numpy.get_include())')" \
    -D BOOST_ROOT="$CONDA_PREFIX" \
    -D Boost_NO_SYSTEM_PATHS=ON \
    -D Boost_NO_BOOST_CMAKE=ON \
    -D RDK_BUILD_AVALON_SUPPORT=ON \
    -D RDK_BUILD_CAIRO_SUPPORT=ON \
    -D RDK_BUILD_CPP_TESTS=OFF \
    -D RDK_BUILD_INCHI_SUPPORT=ON \
    -D RDK_BUILD_FREESASA_SUPPORT=ON \
    -D RDK_BUILD_YAEHMOP_SUPPORT=ON \
    -D RDK_INSTALL_INTREE=OFF \
    -D RDK_INSTALL_STATIC_LIBS=OFF \
    -D RDK_OPTIMIZE_POPCNT=ON \
    ..
make -j 8
make install