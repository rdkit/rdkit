FROM ubuntu:24.04
ARG PYTHON_VERSION=3.12
ARG CODEX_VERSION=latest
ARG DEBIAN_FRONTEND=noninteractive

# System packages: boost 1.83, cairo, eigen, build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    checkinstall \
    wget \
    curl \
    git \
    ca-certificates \
    lsb-release \
    make \
    ninja-build \
    gdb \
    valgrind \
    pkg-config \
    libboost1.83-all-dev \
    libcairo2-dev \
    libeigen3-dev \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (conda with conda-forge as default channel)
RUN ARCH="$(uname -m)" \
    && curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" -o /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:${PATH}"

# Configure conda
RUN conda config --set always_yes yes --set changeps1 no \
    && conda config --set channel_priority strict \
    && conda config --add channels conda-forge

# Create the rdkit_build conda environment (python + python deps only)
RUN conda create --name rdkit_build \
    python=${PYTHON_VERSION} \
    cmake \
    numpy=1.26.4 \
    pillow \
    eigen \
    pandas=2.1 \
    matplotlib-base=3.8 \
    cairo \
    nbval \
    "ipykernel>=6.0" \
    pytest \
    && conda clean -afy

# Activate rdkit_build env via PATH/ENV instead of shell startup
ENV PATH=/opt/conda/envs/rdkit_build/bin:/opt/conda/bin:${PATH}
ENV CONDA_DEFAULT_ENV=rdkit_build
ENV CONDA_PREFIX=/opt/conda/envs/rdkit_build
ENV RDBASE=/workspaces/rdkit-cdd

# Overwrite shell init files to guarantee nothing blocks
# the devcontainer env probe or VS Code server startup.
RUN printf '[ -z "$PS1" ] && return\nsource /opt/conda/etc/profile.d/conda.sh\nconda activate rdkit_build\n' > /root/.bashrc \
    && grep -v conda /etc/bash.bashrc > /tmp/clean.bashrc && mv /tmp/clean.bashrc /etc/bash.bashrc || true \
    && grep -v conda /etc/profile > /tmp/clean.profile && mv /tmp/clean.profile /etc/profile || true \
    && rm -f /etc/profile.d/conda.sh 2>/dev/null || true

# Set up npm global directory and install Codex
ENV NPM_CONFIG_PREFIX=/usr/local/share/npm-global
ENV PATH="${PATH}:/usr/local/share/npm-global/bin"

RUN mkdir -p /usr/local/share/npm-global \
    && npm install -g "@openai/codex@${CODEX_VERSION}" \
    && npm cache clean --force

# Verify the environment
RUN python --version && cmake --version && codex --version
