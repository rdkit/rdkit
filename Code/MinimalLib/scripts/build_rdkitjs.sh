set -e

# Clean and create distribution folder
MINIMALLIB_OUTPUT_PATH="Code/MinimalLib/dist"
rm -rf $MINIMALLIB_OUTPUT_PATH
mkdir -p $MINIMALLIB_OUTPUT_PATH

# Build distribution files
RDKIT_BRANCH=${1:-master}
echo "Building distribution files for release $RDKIT_BRANCH"
DOCKER_BUILDKIT=1 docker build --no-cache -f Code/MinimalLib/docker/Dockerfile --build-arg RDKIT_BRANCH=$RDKIT_BRANCH -o $MINIMALLIB_OUTPUT_PATH .

# Make files executable
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.js
chmod a+rwx $MINIMALLIB_OUTPUT_PATH/RDKit_minimal.wasm

# Log build completed
echo "Build completed"
echo "MinimalLib distribution files are at $MINIMALLIB_OUTPUT_PATH"
