# Example usage of this Dockerfile:
# (the --build-arg arguments are all optional)
#
# 1. cd to Code/MinimalLib/docker
# cd Code/MinimalLib/docker
#
# 2. build the MinimalLib rdkit-minimallib-rdkit-src image:
# docker build --target clone-stage -t rdkit-minimallib-rdkit-src \
#   --build-arg http_proxy=$HTTP_PROXY \
#   --build-arg https_proxy=$HTTP_PROXY \
#   --build-arg "RDKIT_GIT_URL=https://github.com/myfork/rdkit.git" \
#   --build-arg "RDKIT_BRANCH=mybranch" \
#   --network=host -f Dockerfile_rdkit_clone_from_github .


ARG RDKIT_GIT_URL="https://github.com/rdkit/rdkit.git"
ARG RDKIT_BRANCH="master"

FROM rdkit-minimallib-deps AS clone-stage
ARG RDKIT_GIT_URL
ARG RDKIT_BRANCH

LABEL maintainer="Greg Landrum <greg.landrum@t5informatics.com>"

WORKDIR /src
ENV RDBASE=/src/rdkit
RUN git clone -b ${RDKIT_BRANCH} --depth 1 --single-branch ${RDKIT_GIT_URL} 
WORKDIR $RDBASE
