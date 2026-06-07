# Compile WASM MinimalLib

ARG EXCEPTION_HANDLING="-fwasm-exceptions"
ARG FREETYPE_VERSION="2.13.3"

# ---------------------------------------------------------------------------
# Stage 1: build dependencies (emscripten, FreeType, zlib)
# ---------------------------------------------------------------------------
FROM debian:trixie AS deps-stage
ARG EXCEPTION_HANDLING
ARG FREETYPE_VERSION

LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y && apt install -y \
  cmake \
  wget \
  git \
  python3 \
  g++ \
  emscripten \
  libboost-dev \
  libeigen3-dev \
  nodejs

WORKDIR /src
RUN wget -q https://download.savannah.gnu.org/releases/freetype/freetype-${FREETYPE_VERSION}.tar.gz && \
  tar xzf freetype-${FREETYPE_VERSION}.tar.gz
WORKDIR /src/freetype-${FREETYPE_VERSION}
RUN mkdir build
WORKDIR /src/freetype-${FREETYPE_VERSION}/build
RUN emcmake cmake -DCMAKE_BUILD_TYPE=Release \
  -DFT_DISABLE_ZLIB=TRUE -DFT_DISABLE_BZIP2=TRUE -DFT_DISABLE_PNG=TRUE \
  -DFT_DISABLE_HARFBUZZ=TRUE -DFT_DISABLE_BROTLI=TRUE \
  -DCMAKE_C_FLAGS="${EXCEPTION_HANDLING}" -DCMAKE_EXE_LINKER_FLAGS="${EXCEPTION_HANDLING}" \
  -DCMAKE_INSTALL_PREFIX=/opt/freetype ..
RUN make -j2 && make -j2 install

# Pre-build emscripten zlib port into a writable cache
ENV EM_FROZEN_CACHE=0
ENV EM_CACHE=/opt/emscripten-cache
RUN embuilder build zlib

# Copy boost headers and cmake configs into emscripten sysroot,
# patching include paths to avoid host glibc conflicts
RUN cp -r /usr/include/boost ${EM_CACHE}/sysroot/include/ && \
    cp -r /usr/lib/x86_64-linux-gnu/cmake/Boost-* ${EM_CACHE}/sysroot/lib/cmake/ && \
    cp -r /usr/lib/x86_64-linux-gnu/cmake/boost_headers-* ${EM_CACHE}/sysroot/lib/cmake/ && \
    find ${EM_CACHE}/sysroot/lib/cmake/Boost-* ${EM_CACHE}/sysroot/lib/cmake/boost_headers-* \
      -name "*.cmake" -exec sed -i "s|/usr/include|${EM_CACHE}/sysroot/include|g" {} \;

# ---------------------------------------------------------------------------
# Stage 2: copy local rdkit source into the image
# (build context must be the rdkit repo root)
# ---------------------------------------------------------------------------
FROM deps-stage AS src-stage

WORKDIR /
COPY Code /src/rdkit/Code
COPY External /src/rdkit/External
COPY CMakeLists.txt license.txt *.in *.md *.cmake /src/rdkit/

# ---------------------------------------------------------------------------
# Stage 3: configure, patch, build, test
# ---------------------------------------------------------------------------
FROM src-stage AS build-stage
ARG EXCEPTION_HANDLING

WORKDIR /src
ENV RDBASE=/src/rdkit
RUN mkdir build
WORKDIR $RDBASE/build

RUN BOOST_CMAKE_DIR=$(find ${EM_CACHE}/sysroot/lib/cmake -name BoostConfig.cmake | head -1 | xargs dirname) && \
    BOOST_HEADERS_CMAKE_DIR=$(find ${EM_CACHE}/sysroot/lib/cmake -name 'boost_headers-config.cmake' | head -1 | xargs dirname) && \
    emcmake cmake \
  -DRDK_BUILD_FREETYPE_SUPPORT=ON \
  -DRDK_BUILD_MINIMAL_LIB=ON \
  -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
  -DRDK_BUILD_CPP_TESTS=OFF \
  -DRDK_BUILD_INCHI_SUPPORT=ON \
  -DRDK_USE_BOOST_SERIALIZATION=OFF \
  -DRDK_OPTIMIZE_POPCNT=OFF \
  -DRDK_BUILD_THREADSAFE_SSS=OFF \
  -DRDK_BUILD_DESCRIPTORS3D=OFF \
  -DRDK_TEST_MULTITHREADED=OFF \
  -DRDK_BUILD_CHEMDRAW_SUPPORT=OFF \
  -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
  -DRDK_BUILD_COORDGEN_SUPPORT=ON \
  -DRDK_BUILD_MINIMAL_LIB_MCS=ON \
  -DRDK_BUILD_MINIMAL_LIB_MOLZIP=ON \
  -DRDK_BUILD_MINIMAL_LIB_MMPA=ON \
  -DRDK_BUILD_MINIMAL_LIB_RXN=ON \
  -DBoost_DIR="${BOOST_CMAKE_DIR}" \
  -Dboost_headers_DIR="${BOOST_HEADERS_CMAKE_DIR}" \
  -DRDK_BUILD_SLN_SUPPORT=OFF \
  -DRDK_USE_BOOST_IOSTREAMS=OFF \
  -DFREETYPE_INCLUDE_DIRS=/opt/freetype/include/freetype2 \
  -DFREETYPE_LIBRARY=/opt/freetype/lib/libfreetype.a \
  -DZLIB_INCLUDE_DIR="${EM_CACHE}/sysroot/include" \
  -DZLIB_LIBRARY="${EM_CACHE}/sysroot/lib/wasm32-emscripten/libz.a" \
  -DCMAKE_CXX_FLAGS="${EXCEPTION_HANDLING} -O3 -DNDEBUG" \
  -DCMAKE_C_FLAGS="${EXCEPTION_HANDLING} -O3 -DNDEBUG -DCOMPILE_ANSI_ONLY" \
  "-DCMAKE_EXE_LINKER_FLAGS=${EXCEPTION_HANDLING} -s STACK_OVERFLOW_CHECK=1 -s USE_PTHREADS=0 -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s MODULARIZE=1 -s EXPORT_NAME='initRDKitModule' -s USE_ZLIB=1" \
  ..

# Patch to make InChI code work with emscripten
RUN cp /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c \
       /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c.bak && \
  sed 's/&& defined(__APPLE__)//' \
      /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c.bak \
    > /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c

# Patch RapidJSON document.h compilation error (only if the file exists;
# in some rdkit versions rapidjson is fetched later during make)
RUN f=/src/rdkit/External/rapidjson-1.1.0/include/rapidjson/document.h; \
  if [ -f "$f" ]; then \
    sed -i 's|^\( *\)\(GenericStringRef\& operator=(const GenericStringRef\& rhs) { s = rhs.s; length = rhs.length; } *\)$|\1//\2|' "$f"; \
  else \
    echo "rapidjson document.h not found yet, skipping pre-patch (will be patched if needed at build time)"; \
  fi

# Build
RUN make -j2 RDKit_minimal && \
  cp Code/MinimalLib/RDKit_minimal.* ../Code/MinimalLib/demo/

# Run tests
WORKDIR /src/rdkit/Code/MinimalLib/tests
RUN node tests.js

# ---------------------------------------------------------------------------
# Stage 4: export artifacts only (requires BuildKit --output)
# ---------------------------------------------------------------------------
FROM scratch AS export-stage
COPY --from=build-stage /src/rdkit/Code/MinimalLib/demo /
COPY --from=build-stage /src/rdkit/Code/MinimalLib/docs /
