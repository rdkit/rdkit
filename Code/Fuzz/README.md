## Important Notice
The files in the corpora folders (i.e. the folders ending in `_fuzzer`) can not be directly used for purposes other than fuzzing.
This is because the fuzzer uses parts of the content for generating different information.
Consider the example `[OH3+]0`.
The first part `[OH3+]` will be used as a smiles formula, but the last part `0` will for example be used to determine
whether the fuzzer should set a certain flag to `true` or it will be used to derive an integral value.

## Compiling
To fuzz rdkit you need to have clang installed.
If you have built the fuzzers you can invoke them like this:
./smiles_string_to_mol_fuzzer -dict=smiles_string_to_mol_fuzzer.dict smiles_string_to_mol_fuzzer/
For possible options that you can pass to the fuzzer see the libFuzzer [docs](https://llvm.org/docs/LibFuzzer.html).
# Clang
````shell
export CC="clang"
export CXX="clang++"
export SANITIZER_FLAGS_address="-fsanitize=address -fsanitize-address-use-after-scope"
export COVERAGE_FLAGS="-fsanitize=fuzzer-no-link"
export CFLAGS="-O1 -fno-omit-frame-pointer -gline-tables-only -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION $COVERAGE_FLAGS $SANITIZER_FLAGS_address"
export CXXFLAGS="$CFLAGS"
export LIB_FUZZING_ENGINE="-fsanitize=fuzzer"

mkdir build && cd build && \
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DRDK_INSTALL_INTREE=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DLIB_FUZZING_ENGINE=${LIB_FUZZING_ENGINE} \
    -DRDK_BUILD_FUZZ_TARGETS=ON \
    -DRDK_INSTALL_STATIC_LIBS=ON \
    -DBoost_USE_STATIC_LIBS=ON \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DBoost_NO_SYSTEM_PATHS=ON \
make
````

# GCC (non-fuzzing mode)
In this mode the resulting fuzzers take a list of files as argument
and invoke the fuzz target on each file.
No actual fuzzing will be done, since no new test cases are generated.
````shell
export CC="gcc"
export CXX="g++"
export SANITIZER_FLAGS_address="-fsanitize=address -fsanitize-address-use-after-scope"
export COVERAGE_FLAGS=""
export CFLAGS="-O1 -fno-omit-frame-pointer -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION $COVERAGE_FLAGS $SANITIZER_FLAGS_address"
export CXXFLAGS="$CFLAGS"
export LIB_FUZZING_ENGINE=""

mkdir build && cd build && \
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DRDK_INSTALL_INTREE=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DLIB_FUZZING_ENGINE=${LIB_FUZZING_ENGINE} \
    -DRDK_BUILD_FUZZ_TARGETS=ON \
    -DRDK_INSTALL_STATIC_LIBS=ON \
    -DBoost_USE_STATIC_LIBS=ON \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DBoost_NO_SYSTEM_PATHS=ON \
make
````

# GCC (fuzzing mode)
This does not seem to be possible.