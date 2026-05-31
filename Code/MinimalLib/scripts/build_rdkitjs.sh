#!/usr/bin/env bash

set -e

# Resolve MinimalLib root (two levels up from this script)
cwd=$(realpath $0)
cwd=$(dirname $cwd)
cwd=$(dirname $cwd)
MINIMALLIB_ROOT=$cwd

# Resolve repo root (two levels up from MinimalLib: Code/MinimalLib -> Code -> repo root)
REPO_ROOT=$(dirname $(dirname $MINIMALLIB_ROOT))

# Clean and create distribution folder
MINIMALLIB_OUTPUT_PATH=$REPO_ROOT/dist
rm -rf $MINIMALLIB_OUTPUT_PATH
mkdir -p $MINIMALLIB_OUTPUT_PATH

# Build distribution files
if [ "$GET_SRC" = copy_from_local ]; then
    echo "Building distribution files from local source tree"
    docker buildx bake \
        -f $MINIMALLIB_ROOT/docker/docker-bake.hcl \
        --set "export.output=type=local,dest=$MINIMALLIB_OUTPUT_PATH" \
        export
else
    RDKIT_BRANCH=${1:-master}
    RDKIT_GIT_URL=${2:-"https://github.com/rdkit/rdkit.git"}
    echo "Building distribution files for release $RDKIT_BRANCH from repo $RDKIT_GIT_URL"
    docker compose -f $MINIMALLIB_ROOT/docker/docker_compose_build_minimallib.yml build \
        --no-cache \
        --build-arg "EXCEPTION_HANDLING=-fwasm-exceptions" \
        --build-arg "RDKIT_GIT_URL=$RDKIT_GIT_URL" \
        --build-arg "RDKIT_GIT_BRANCH=$RDKIT_BRANCH"
    DOCKER_BUILDKIT=1 docker build \
        -f $MINIMALLIB_ROOT/docker/Dockerfile_4_rdkit_export \
        -o $MINIMALLIB_OUTPUT_PATH \
        $MINIMALLIB_ROOT/docker
fi

# Make files executable
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.js
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.wasm

# Log build completed
echo "Build completed"
echo "MinimalLib distribution files are at $MINIMALLIB_OUTPUT_PATH"
