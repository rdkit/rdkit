#!/usr/bin/env bash

set -e

GET_SRC=${GET_SRC:-copy}
NPM_VERSION=${NPM_VERSION:-0.0.0}

# Set buildkit to false to force sequential docker compilation
export DOCKER_BUILDKIT=0

# Clean and create distribution folder
cwd=$(realpath $0)
cwd=$(dirname $cwd)
cwd=$(dirname $cwd)
cd $cwd
MINIMALLIB_OUTPUT_PATH=$(realpath build)
rm -rf $MINIMALLIB_OUTPUT_PATH
mkdir -p $MINIMALLIB_OUTPUT_PATH

# Build distribution files
if [ "$GET_SRC" = copy ]; then
    echo "Building distribution files from local source tree"
    docker compose -f docker/docker_compose_build_minimallib.yml build \
        --build-arg "EXCEPTION_HANDLING=-fwasm-exceptions" \
        --build-arg "VERSION=${NPM_VERSION}"
else
    RDKIT_BRANCH=${1:-master}
    RDKIT_GIT_URL=${2:-"https://github.com/rdkit/rdkit.git"}
    echo "Building distribution files for release $RDKIT_BRANCH from repo $RDKIT_GIT_URL"
    GET_SRC=clone docker compose -f docker/docker_compose_build_minimallib.yml build \
        --no-cache \
        --build-arg "EXCEPTION_HANDLING=-fwasm-exceptions" \
        --build-arg "RDKIT_GIT_URL=$RDKIT_GIT_URL" \
        --build-arg "RDKIT_GIT_BRANCH=$RDKIT_BRANCH"
fi
DOCKER_BUILDKIT=1 docker build -f docker/04_export.Dockerfile -o $MINIMALLIB_OUTPUT_PATH .

echo "Build completed"
echo "MinimalLib distribution files are at $MINIMALLIB_OUTPUT_PATH"
