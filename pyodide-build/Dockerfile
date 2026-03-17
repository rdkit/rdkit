##############################################################################
# Dockerfile: Build RDKit Python wrappers for Pyodide (wasm32-emscripten)
#
# Build (from rdkit repo root):
#   docker build -t rdkit-pyodide -f pyodide-build/Dockerfile .
#
# Extract wheel:
#   docker create --name tmp rdkit-pyodide
#   docker cp tmp:/out/rdkit_pyodide-*.whl ./pyodide-build/
#   docker rm tmp
#
# Test (requires node and npm):
#   cd pyodide-build
#   npm install
#   node test_pyodide.js
##############################################################################

FROM python:3.13-bookworm

# Pyodide 0.29.3: Emscripten 4.0.9, Python 3.13, wasm-exceptions
ARG EMSDK_VERSION="4.0.9"
ARG BOOST_VERSION="1.85.0"
ARG BOOST_VERSION_U="1_85_0"
ARG ZLIB_VERSION="1.3.2"
ARG PYODIDE_BUILD_VERSION="0.29.3"
ARG PYTHON_VERSION="3.13"
# Pyodide 0.29+ uses native WASM exceptions
ARG EXCEPTION_HANDLING="-fwasm-exceptions"

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies + numpy (for C API headers needed by boost.numpy)
RUN apt-get update && apt-get install -y \
    curl wget g++ git xz-utils libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install --break-system-packages "cmake==3.31.6" "numpy>=2"

SHELL ["/bin/bash", "-c"]

# Get numpy include path (host numpy, just for headers)
RUN python3 -c "import numpy; print(numpy.get_include())" > /tmp/numpy_include.txt

# ========================
# 1. Install Emscripten SDK
# ========================
WORKDIR /opt
RUN git clone https://github.com/emscripten-core/emsdk.git
WORKDIR /opt/emsdk
RUN ./emsdk install ${EMSDK_VERSION} && ./emsdk activate ${EMSDK_VERSION}
RUN echo "source /opt/emsdk/emsdk_env.sh > /dev/null 2>&1" >> ~/.bashrc
ENV EMSDK=/opt/emsdk
ENV PATH="/opt/emsdk:/opt/emsdk/upstream/emscripten:${PATH}"

# ========================
# 2. Install pyodide-build and get xbuildenv (for CPython headers)
# ========================
RUN pip3 install --break-system-packages "pyodide-build==${PYODIDE_BUILD_VERSION}" "wheel>=0.40,<0.46"

# Download xbuildenv to get cross-compiled CPython headers
RUN pyodide xbuildenv install --path /opt/xbuildenv ${PYODIDE_BUILD_VERSION}

# Find the Python include directory dynamically
RUN find /opt/xbuildenv -name "Python.h" -exec dirname {} \; > /tmp/python_include.txt && \
    cat /tmp/python_include.txt

# ========================
# 3. Build zlib for Emscripten
# ========================
WORKDIR /src
RUN wget -q "https://github.com/madler/zlib/releases/download/v${ZLIB_VERSION}/zlib-${ZLIB_VERSION}.tar.gz" \
    && tar xf "zlib-${ZLIB_VERSION}.tar.gz"
WORKDIR /src/zlib-${ZLIB_VERSION}
RUN source /opt/emsdk/emsdk_env.sh && \
    CFLAGS="${EXCEPTION_HANDLING}" \
    emconfigure ./configure --prefix=/opt/zlib --static && \
    emmake make -j$(nproc) && emmake make install

# ========================
# 4. Build Boost for Emscripten (with boost.python + boost.numpy!)
# ========================
WORKDIR /src
RUN wget -q "https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION_U}.tar.bz2" \
    && tar xf "boost_${BOOST_VERSION_U}.tar.bz2"
WORKDIR /src/boost_${BOOST_VERSION_U}

# Patch boost.numpy for NumPy 2.x compatibility (elsize -> PyDataType_ELSIZE)
RUN sed -i 's/reinterpret_cast<PyArray_Descr\*>(ptr())->elsize/PyDataType_ELSIZE(reinterpret_cast<PyArray_Descr*>(ptr()))/' \
    libs/python/src/numpy/dtype.cpp

# Bootstrap boost with host toolchain
RUN ./bootstrap.sh --prefix=/opt/boost --with-libraries=system,python,serialization,iostreams

# Create user-config.jam for emscripten cross-compilation
# Include both Python and NumPy headers so boost.numpy can build
RUN PYINC=$(cat /tmp/python_include.txt) && \
    NUMPY_INCLUDE=$(cat /tmp/numpy_include.txt) && \
    PYVER=$(echo ${PYTHON_VERSION} | cut -d. -f1-2) && \
    printf '%s\n' \
    "using gcc : emscripten : /opt/emsdk/upstream/emscripten/em++ :" \
    "    <cxxflags>\"${EXCEPTION_HANDLING} -std=c++20 -fPIC\"" \
    "    <linkflags>\"${EXCEPTION_HANDLING}\"" \
    "    <archiver>/opt/emsdk/upstream/emscripten/emar" \
    "    <ranlib>/opt/emsdk/upstream/emscripten/emranlib" \
    "    ;" \
    "using python : ${PYVER} : : ${PYINC} : ${NUMPY_INCLUDE} ;" \
    > user-config.jam && cat user-config.jam

# Build boost — static only, single-threaded, 32-bit (wasm32)
RUN source /opt/emsdk/emsdk_env.sh && \
    ./b2 install \
    --prefix=/opt/boost \
    --user-config=user-config.jam \
    toolset=gcc-emscripten \
    variant=release \
    link=static \
    threading=single \
    runtime-link=static \
    address-model=32 \
    --with-system \
    --with-python \
    --with-serialization \
    --with-iostreams \
    -sZLIB_INCLUDE=/opt/zlib/include \
    -sZLIB_LIBPATH=/opt/zlib/lib \
    -j$(nproc) 2>&1 ; true

# Verify boost.python and boost.numpy were built
RUN echo "=== Boost libraries ===" && ls -la /opt/boost/lib/*.a 2>/dev/null && \
    echo "=== Boost cmake configs ===" && ls /opt/boost/lib/cmake/

# ========================
# 5. Build RDKit — all STATIC (no shared linking issues with Emscripten)
# ========================
WORKDIR /src
# Copy local RDKit source tree (build context is rdkit repo root)
COPY . /src/rdkit

ENV RDBASE=/src/rdkit

WORKDIR /src/rdkit/build

# Create a dummy libpython archive — Emscripten side modules don't
# actually link against libpython (it's part of the main Pyodide module)
# but CMake's FindPython3 needs it to create the Python3::Python target.
RUN source /opt/emsdk/emsdk_env.sh && \
    PYINC=$(cat /tmp/python_include.txt) && \
    PYVER=$(echo ${PYTHON_VERSION} | cut -d. -f1-2) && \
    mkdir -p /opt/python-dummy/lib && \
    emar rcs /opt/python-dummy/lib/libpython${PYVER}.a && \
    NUMPY_INCLUDE=$(cat /tmp/numpy_include.txt) && \
    emcmake cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DRDK_BUILD_PYTHON_WRAPPERS=ON \
    -DRDK_BUILD_MINIMAL_LIB=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_USE_BOOST_SERIALIZATION=ON \
    -DRDK_USE_BOOST_IOSTREAMS=ON \
    -DRDK_OPTIMIZE_POPCNT=OFF \
    -DRDK_BUILD_THREADSAFE_SSS=OFF \
    -DRDK_BUILD_DESCRIPTORS3D=ON \
    -DRDK_TEST_MULTITHREADED=OFF \
    -DRDK_BUILD_CHEMDRAW_SUPPORT=OFF \
    -DRDK_BUILD_MAEPARSER_SUPPORT=OFF \
    -DRDK_BUILD_COORDGEN_SUPPORT=OFF \
    -DRDK_BUILD_SLN_SUPPORT=OFF \
    -DRDK_BUILD_CAIRO_SUPPORT=OFF \
    -DRDK_BUILD_QT_SUPPORT=OFF \
    -DRDK_BUILD_PGSQL=OFF \
    -DRDK_BUILD_SWIG_WRAPPERS=OFF \
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_BUILD_AVALON_SUPPORT=OFF \
    -DRDK_BUILD_YAEHMOP_SUPPORT=OFF \
    -DRDK_BUILD_XYZ2MOL_SUPPORT=OFF \
    -DRDK_INSTALL_INTREE=ON \
    -DRDK_INSTALL_STATIC_LIBS=ON \
    -DBoost_USE_STATIC_LIBS=ON \
    -DBoost_USE_STATIC_RUNTIME=ON \
    -DCMAKE_PREFIX_PATH="/opt/boost;/opt/zlib" \
    -DCMAKE_FIND_ROOT_PATH="/opt/boost;/opt/zlib;/opt/xbuildenv" \
    -DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH \
    -DBoost_DIR=/opt/boost/lib/cmake/Boost-${BOOST_VERSION} \
    -Dboost_headers_DIR=/opt/boost/lib/cmake/boost_headers-${BOOST_VERSION} \
    -Dboost_python_DIR=/opt/boost/lib/cmake/boost_python-${BOOST_VERSION} \
    -Dboost_numpy_DIR=/opt/boost/lib/cmake/boost_numpy-${BOOST_VERSION} \
    -DBoost_INCLUDE_DIR=/opt/boost/include \
    -DZLIB_INCLUDE_DIR=/opt/zlib/include \
    -DZLIB_LIBRARY=/opt/zlib/lib/libz.a \
    -DPython3_INCLUDE_DIR=${PYINC} \
    -DPython3_LIBRARY=/opt/python-dummy/lib/libpython${PYVER}.a \
    -DPython3_NumPy_INCLUDE_DIR=${NUMPY_INCLUDE} \
    -DPython3_EXECUTABLE=/usr/local/bin/python3 \
    -DCMAKE_CXX_FLAGS="${EXCEPTION_HANDLING} -O2 -DNDEBUG -fPIC -std=c++20" \
    -DCMAKE_C_FLAGS="${EXCEPTION_HANDLING} -O2 -DNDEBUG -fPIC" \
    .. 2>&1 | tee /tmp/cmake_output.txt

# Build everything (all static — MODULE targets become .a due to Emscripten)
RUN source /opt/emsdk/emsdk_env.sh && \
    make -j$(nproc) 2>&1 | tee /tmp/build_output.txt

# ========================
# 6. Relink Python wrapper .a archives into .so side modules
# ========================
# The cmake build produces static .a archives for the Python extension modules
# because Emscripten disables shared libs by default. We relink each wrapper
# into a .so WASM side module.
#
# Strategy: Create one "core" .so with all RDKit C++ libraries (loaded first),
# then small wrapper .so modules that import symbols from the core at runtime.
# Pyodide's dlopen adds exported symbols to the global namespace, so later
# modules can resolve against earlier ones.
RUN source /opt/emsdk/emsdk_env.sh && \
    mkdir -p /tmp/relink && \
    echo "=== Collecting all static libraries ===" && \
    # Collect all .a files (RDKit + Boost + zlib)
    find /src/rdkit/build/lib -name '*.a' | sort > /tmp/relink/all_libs.txt && \
    echo "/opt/boost/lib/libboost_python313.a" >> /tmp/relink/all_libs.txt && \
    echo "/opt/boost/lib/libboost_numpy313.a" >> /tmp/relink/all_libs.txt && \
    echo "/opt/boost/lib/libboost_system.a" >> /tmp/relink/all_libs.txt && \
    echo "/opt/boost/lib/libboost_serialization.a" >> /tmp/relink/all_libs.txt && \
    echo "/opt/boost/lib/libboost_iostreams.a" >> /tmp/relink/all_libs.txt && \
    echo "/opt/zlib/lib/libz.a" >> /tmp/relink/all_libs.txt && \
    echo "Total libraries: $(wc -l < /tmp/relink/all_libs.txt)" && \
    echo "=== Building core side module ===" && \
    # Link all archives with --whole-archive to include every object file,
    # --export-all to export every symbol, --no-gc-sections to prevent stripping
    WHOLE_ARGS="" && \
    while read -r lib; do \
        WHOLE_ARGS="$WHOLE_ARGS -Wl,--whole-archive $lib -Wl,--no-whole-archive"; \
    done < /tmp/relink/all_libs.txt && \
    em++ -fwasm-exceptions -sSIDE_MODULE=1 -O2 -shared \
        -Wl,--no-gc-sections -Wl,--export-all \
        $WHOLE_ARGS \
        -o /src/rdkit/build/librdkit_core.so 2>&1 && \
    ls -lh /src/rdkit/build/librdkit_core.so && \
    echo "=== Relinking Python wrappers as .so side modules ===" && \
    # Each wrapper .so only contains its own code; symbols resolve from core
    for wrapper_a in $(find /src/rdkit/build -path '*/Wrap/*' -name '*.a' ! -name '*_static*'); do \
        modname=$(basename "$wrapper_a" .a) && \
        wrapper_dir=$(dirname "$wrapper_a") && \
        obj_files=$(find "$wrapper_dir/CMakeFiles/${modname}.dir" -name '*.o' 2>/dev/null) && \
        if [ -z "$obj_files" ]; then \
            echo "SKIP: $modname (no .o files found)" && \
            continue; \
        fi && \
        out_so="/src/rdkit/build/${modname}.so" && \
        echo "LINK: $modname" && \
        em++ -fwasm-exceptions -sSIDE_MODULE=1 -O2 -shared \
            $obj_files \
            -o "$out_so" 2>&1 || echo "FAILED: $modname"; \
    done && \
    echo "=== Side module relinking complete ===" && \
    find /src/rdkit/build -name "*.so" | wc -l && echo "total .so files"

# ========================
# 7. Collect outputs into proper package structure
# ========================
# Copy Python source files first, then place .so files in correct subdirs
RUN mkdir -p /out/rdkit && \
    cp -r /src/rdkit/rdkit/* /out/rdkit/ 2>/dev/null || true && \
    # Copy core module
    cp /src/rdkit/build/librdkit_core.so /out/rdkit/ && \
    # Place each wrapper .so in its correct package directory
    # (mapping from rdkit_python_extension DEST parameters)
    for so_file in /src/rdkit/build/*.so; do \
        modname=$(basename "$so_file" .so) && \
        case "$modname" in \
            librdkit_core) ;; \
            rdBase) cp "$so_file" /out/rdkit/ ;; \
            cDataStructs) cp "$so_file" /out/rdkit/DataStructs/ ;; \
            DistGeom) cp "$so_file" /out/rdkit/DistanceGeometry/ ;; \
            rdForceField) cp "$so_file" /out/rdkit/ForceField/ ;; \
            rdGeometry) cp "$so_file" /out/rdkit/Geometry/ ;; \
            rdSimDivPickers) cp "$so_file" /out/rdkit/SimDivFilters/ ;; \
            rdMetricMatrixCalc) mkdir -p /out/rdkit/DataManip/Metric && cp "$so_file" /out/rdkit/DataManip/Metric/ ;; \
            rdInfoTheory) mkdir -p /out/rdkit/ML/InfoTheory && cp "$so_file" /out/rdkit/ML/InfoTheory/ ;; \
            rdAlignment) cp "$so_file" /out/rdkit/Numerics/ ;; \
            rdMolDraw2D) cp "$so_file" /out/rdkit/Chem/Draw/ ;; \
            rdMolStandardize) mkdir -p /out/rdkit/Chem/MolStandardize && cp "$so_file" /out/rdkit/Chem/MolStandardize/ ;; \
            rdScaffoldNetwork) mkdir -p /out/rdkit/Chem/Scaffolds && cp "$so_file" /out/rdkit/Chem/Scaffolds/ ;; \
            *) cp "$so_file" /out/rdkit/Chem/ ;; \
        esac; \
    done && \
    # Copy data files inside the rdkit package so they're always available
    cp -r /src/rdkit/Data /out/rdkit/Data 2>/dev/null || true && \
    # Generate RDPaths.py that resolves paths relative to the package directory
    printf '%s\n' \
        'import os' \
        '_share = os.path.dirname(__file__)' \
        'RDCodeDir = _share' \
        'RDDataDir = os.path.join(_share, "Data")' \
        'RDDocsDir = os.path.join(_share, "Docs")' \
        'RDProjDir = os.path.join(_share, "Projects")' \
        'RDContribDir = os.path.join(_share, "Contrib")' \
        > /out/rdkit/RDPaths.py && \
    echo "=== Renaming .so to .so.wasm to prevent micropip auto-loading ===" && \
    find /out/rdkit -name "*.so" -exec sh -c 'mv "$1" "$1.wasm"' _ {} \; && \
    echo "=== __init__.py will be patched in next step ===" && \
    echo "=== Build complete ===" && \
    find /out/rdkit -name "*.so.wasm" | sort && \
    echo "---" && \
    du -sh /out/rdkit

# ========================
# 8. Package as Pyodide wheel
# ========================
RUN cat > /tmp/patch_init.py << 'PYEOF'
import sys
init = open('/out/rdkit/__init__.py').read()
loader = '''import sys as _sys

if _sys.platform == 'emscripten':
    import os as _os
    import importlib.abc
    import importlib.machinery

    from pyodide_js._module import loadDynamicLibrary as _ldl
    import js as _js
    _core = _os.path.join(_os.path.dirname(__file__), 'librdkit_core.so.wasm')
    _ldl(_core, _js.JSON.parse('{"global": true}'))
    del _ldl, _js, _core

    class _RDKitExtensionFinder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path, target=None):
            parts = fullname.split('.')
            if parts[0] != 'rdkit':
                return None
            modname = parts[-1]
            if path:
                for d in path:
                    candidate = _os.path.join(d, modname + '.so.wasm')
                    if _os.path.exists(candidate):
                        loader = importlib.machinery.ExtensionFileLoader(
                            fullname, candidate
                        )
                        return importlib.util.spec_from_file_location(
                            fullname, candidate, loader=loader,
                        )
            return None

    import importlib.util
    _sys.meta_path.insert(0, _RDKitExtensionFinder())

'''
open('/out/rdkit/__init__.py', 'w').write(loader + init)
print('Patched __init__.py')
PYEOF
RUN python3 /tmp/patch_init.py

RUN cat > /tmp/build_wheel.py << 'PYEOF'
import os, re, zipfile
with open("/src/rdkit/CMakeLists.txt") as f:
    txt = f.read()
m = re.search(r'set\(RDKit_Year "(\d+)"\).*set\(RDKit_Month "(\d+)"\).*set\(RDKit_Revision "(\d+)"\)', txt, re.S)
if not m:
    raise RuntimeError("Could not parse RDKit version")
version = f"{m.group(1)}.{m.group(2)}.{m.group(3)}"
print(f"RDKit version: {version}")
whl_tag = "cp313-cp313-pyodide_2025_0_wasm32"
whl_name = f"rdkit_pyodide-{version}-{whl_tag}"
dist_info = f"/out/{whl_name}.dist-info"
os.makedirs(dist_info, exist_ok=True)
open(f"{dist_info}/METADATA", "w").write(f"Metadata-Version: 2.1\nName: rdkit-pyodide\nVersion: {version}\nSummary: RDKit cheminformatics library for Pyodide\n")
open(f"{dist_info}/WHEEL", "w").write(f"Wheel-Version: 1.0\nGenerator: rdkit-pyodide-build\nRoot-Is-Purelib: false\nTag: {whl_tag}\n")
open(f"{dist_info}/top_level.txt", "w").write("rdkit\n")
open(f"{dist_info}/RECORD", "w").close()
whl_path = f"/out/{whl_name}.whl"
with zipfile.ZipFile(whl_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for top_dir in ["rdkit", f"{whl_name}.dist-info"]:
        base = f"/out/{top_dir}"
        for root, _, files in os.walk(base):
            for fname in files:
                full = os.path.join(root, fname)
                zf.write(full, os.path.relpath(full, "/out"))
print(f"Wheel: {whl_path} ({os.path.getsize(whl_path)/1024/1024:.1f} MB)")
PYEOF
RUN python3 /tmp/build_wheel.py && ls -lh /out/*.whl

CMD ["echo", "RDKit Pyodide build complete. Use docker cp to extract /out/*.whl"]
