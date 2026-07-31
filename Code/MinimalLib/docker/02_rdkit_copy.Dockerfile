# Example usage of this Dockerfile:
# (the --build-arg arguments are all optional)
#
# 1. cd to Code/MinimalLib/docker
# cd Code/MinimalLib/docker
#
# 2. build the MinimalLib rdkit-minimallib-rdkit-src image:
# docker build --target local-src-stage -t rdkit-minimallib-rdkit-src \
#   --build-arg http_proxy=$HTTP_PROXY \
#   --build-arg https_proxy=$HTTP_PROXY \
#   --network=host -f Dockerfile_2_rdkit_copy_from_local ../../..


FROM rdkit-minimallib-deps AS local-src-stage

LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

WORKDIR /
COPY Code /src/rdkit/Code
COPY External /src/rdkit/External
COPY CMakeLists.txt license.txt *.in *.md *.cmake /src/rdkit/
