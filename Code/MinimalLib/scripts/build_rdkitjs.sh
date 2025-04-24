set -e

# Clean and create distribution folder
MINIMALLIB_OUTPUT_PATH="Code/MinimalLib/dist"
rm -rf $MINIMALLIB_OUTPUT_PATH
mkdir -p $MINIMALLIB_OUTPUT_PATH

# Build distribution files
RDKIT_BRANCH=${1:-master}
RDKIT_GIT_URL=${2:-"https://github.com/rdkit/rdkit.git"}
echo "Building distribution files for release $RDKIT_BRANCH from repo $RDKIT_GIT_URL"
GET_SRC=clone_from_github docker-compose -f docker/docker_compose_build_minimallib.yml build --no-cache \
  --build-arg "EXCEPTION_HANDLING=-fwasm-exceptions" \
  --build-arg "RDKIT_GIT_URL=$RDKIT_GIT_URL" \
  --build-arg "RDKIT_GIT_BRANCH=$RDKIT_BRANCH"
DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile_4_rdkit_export  -o $MINIMALLIB_OUTPUT_PATH .

# Make files executable
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.js
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.wasm

# Log build completed
echo "Build completed"
echo "MinimalLib distribution files are at $MINIMALLIB_OUTPUT_PATH"
