# Example usage of this Dockerfile:
# (the build-arg arguments are all optional)
#
# 1. cd to Code/MinimalLib/docker
# cd Code/MinimalLib/docker
#
# 2. build the JS and WASM libraries
#    (the build-arg arguments are all optional)
# docker build -t rdkit-minimallib --network=host \
#  --build-arg "RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS=ON" \
#  --build-arg "RDKIT_GIT_URL=https://github.com/myfork/rdkit.git" \
#  --build-arg "RDKIT_BRANCH=mybranch"
#
# 3. create a temporary container and copy built libraries
#    from the container to your local filesystem, then destroy
#    the temporary container
# docker create --name=rdkit-minimallib-container rdkit-minimallib:latest --entrypoint /
# docker cp rdkit-minimallib-container:/RDKit_minimal.js ../demo
# docker cp rdkit-minimallib-container:/RDKit_minimal.wasm ../demo
# docker cp rdkit-minimallib-container:/RDKit_minimal_legacy.js ../demo
# docker rm rdkit-minimallib-container


ARG RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS="OFF"
ARG RDKIT_GIT_URL="https://github.com/rdkit/rdkit.git"
ARG RDKIT_BRANCH="master"

FROM debian:buster as build-stage
ARG RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS
ARG RDKIT_GIT_URL
ARG RDKIT_BRANCH

LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

RUN apt-get update && apt-get upgrade -y && apt install -y \
  ca-certificates \
  curl \
  wget \
  cmake \
  python3 \
  g++ \
  libeigen3-dev \
  git \
  nodejs

ENV LANG C

WORKDIR /opt
RUN wget -q https://boostorg.jfrog.io/artifactory/main/release/1.67.0/source/boost_1_67_0.tar.gz && \
  tar xzf boost_1_67_0.tar.gz 
WORKDIR /opt/boost_1_67_0
RUN ./bootstrap.sh --prefix=/opt/boost --with-libraries=system && \
  ./b2 install


WORKDIR /opt
RUN git clone https://github.com/emscripten-core/emsdk.git

WORKDIR /opt/emsdk
RUN ./emsdk update-tags && \
  ./emsdk install latest && \
  ./emsdk activate latest

#RUN source ./emsdk_env.sh

RUN mkdir /src
WORKDIR /src
ENV RDBASE=/src/rdkit
RUN git clone ${RDKIT_GIT_URL}
WORKDIR $RDBASE
RUN git fetch --all --tags && \
  git checkout ${RDKIT_BRANCH}

# this is to display the emsdk_env.sh output only once per build
RUN echo "[ -e ~/emsdk_env.log ] && source /opt/emsdk/emsdk_env.sh > /dev/null 2>&1 || source /opt/emsdk/emsdk_env.sh 2>&1 | tee ~/emsdk_env.log" >> ~/.bashrc
SHELL ["/bin/bash", "-c", "-l"]
RUN echo -e '#!/usr/bin/env bash\n\n\
emcmake cmake -DBoost_INCLUDE_DIR=/opt/boost/include -DRDK_BUILD_FREETYPE_SUPPORT=ON -DRDK_BUILD_MINIMAL_LIB=ON \
-DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_INCHI_SUPPORT=ON \
-DRDK_USE_BOOST_SERIALIZATION=OFF -DRDK_OPTIMIZE_POPCNT=OFF -DRDK_BUILD_THREADSAFE_SSS=OFF \
-DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_TEST_MULTITHREADED=OFF \
-DRDK_BUILD_MAEPARSER_SUPPORT=OFF -DRDK_BUILD_COORDGEN_SUPPORT=ON \
-DRDK_BUILD_SLN_SUPPORT=OFF -DRDK_USE_BOOST_IOSTREAMS=OFF \
-DFREETYPE_INCLUDE_DIRS=/opt/emsdk/upstream/emscripten/cache/sysroot/include/freetype2 \
-DFREETYPE_LIBRARY=/opt/emsdk/upstream/emscripten/cache/sysroot/lib/wasm32-emscripten/libfreetype.a \
-DCMAKE_CXX_FLAGS="-Wno-enum-constexpr-conversion -s DISABLE_EXCEPTION_CATCHING=0" \
-DCMAKE_C_FLAGS="-Wno-enum-constexpr-conversion -DCOMPILE_ANSI_ONLY" \
-DCMAKE_EXE_LINKER_FLAGS="$EXTRA_LINKER_FLAGS -s MODULARIZE=1 -s EXPORT_NAME=\\"'"'"initRDKitModule"'"'\\"" $@' > /src/emcmake.sh
run chmod +x /src/emcmake.sh

RUN mkdir build
WORKDIR build

# run emcmake
RUN /src/emcmake.sh ..

# "patch" to make the InChI code work with emscripten:
RUN cp /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c.bak && \
  sed 's/&& defined(__APPLE__)//' /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c.bak > /src/rdkit/External/INCHI-API/src/INCHI_BASE/src/util.c

# build and "install"
RUN make -j2 RDKit_minimal && cp Code/MinimalLib/RDKit_minimal.* ../Code/MinimalLib/demo/

# also build the legacy, plain JS version if required
# to make it work on legacy browser we need a few polyfills
# and also to change an arrow function into a standard function
RUN [ "$RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS" = "ON" ] && \
  rm -f CMakeCache.txt && \
    EXTRA_LINKER_FLAGS="-s SINGLE_FILE=1 -s LEGACY_VM_SUPPORT=1 -s MIN_IE_VERSION=11 -s WASM=0 --memory-init-file 0" \
      /src/emcmake.sh -DRDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS=ON .. && \
    make -j2 RDKit_minimal && echo -e '\
var toObj=function(c){if(null===c||"undefined"===typeof c)throw new TypeError("this is null or undefined");return Object(c)};Uint8Array.prototype.fill||Object.defineProperty(Uint8Array.prototype,"fill",{value:function(c,b,d){var e=toObj(this),a=e.length>>>0;b>>=0;b=0>b?Math.max(a+b,0):Math.min(b,a);d="undefined"===typeof d?a:d>>0;for(a=0>d?Math.max(a+d,0):Math.min(d,a);b<a;)e[b]=c,++b;return e}});\n\
Uint8Array.prototype.copyWithin||Object.defineProperty(Uint8Array.prototype,"copyWithin",{value:function(c,b,d){var e=toObj(this),a=e.length>>>0;c>>=0;c=0>c?Math.max(a+c,0):Math.min(c,a);b>>=0;b=0>b?Math.max(a+b,0):Math.min(b,a);d="undefined"===typeof d?a:d>>0;a=Math.min((0>d?Math.max(a+d,0):Math.min(d,a))-b,a-c);d=1;b<c&&c<b+a&&(d=-1,b+=a-1,c+=a-1);for(;0<a;)b in e?e[c]=e[b]:delete e[c],b+=d,c+=d,--a;return e}});\n\
String.prototype.includes||Object.defineProperty(String.prototype,"includes",{value:function(c,b){var d=toObj(this);if(c instanceof RegExp)throw new TypeError("first argument must not be a RegExp");return-1!==d.indexOf(c,b||0)}});\n\
Array.prototype.includes||Object.defineProperty(Array.prototype,"includes",{value:function(c,b){var d=toObj(this),e=d.length>>>0;if(0===e)return!1;var a=b||0;for(a=Math.max(0<=a?a:e-Math.abs(a),0);a<e;){var f=d[a],g=c;if(f===g||"number"===typeof f&&"number"===typeof g&&isNaN(f)&&isNaN(g))return!0;++a}return!1}});String.prototype.startsWith||Object.defineProperty(String.prototype,"startsWith",{value:function(c,b){var d=toObj(this),e=0<b?b|0:0;return d.substring(e,e+c.length)===c}});\n\
' > ../Code/MinimalLib/demo/RDKit_minimal_legacy.js && \
    sed 's|\(var *initRDKitModule *= *(\)\(() *=>\)\( *{\)|\1function()\3|' \
      < Code/MinimalLib/RDKit_minimal.js >> ../Code/MinimalLib/demo/RDKit_minimal_legacy.js

WORKDIR /src/rdkit/Code/MinimalLib/tests
RUN node -v
# run the tests
RUN RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS=OFF node tests.js
# run the tests again on the legacy version if we built it
RUN [ "$RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS" = ON ] && node tests.js

# Copy js and wasm rdkit files to use in browser
# This feature requires the BuildKit backend
# https://docs.docker.com/engine/reference/commandline/build/#custom-build-outputs
FROM scratch as export-stage
COPY --from=build-stage /src/rdkit/Code/MinimalLib/demo /
COPY --from=build-stage /src/rdkit/Code/MinimalLib/docs /
