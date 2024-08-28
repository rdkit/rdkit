# Wrapper to make using 

## Compilation and installation

In order to build the code in this repo, you need to also have the pubchem-align3d code base, which is brought in via a git submodule. So you either need to provide the `--recurse-submodules` argument when you clone this repo or clone it normally and then run:
```
% git submodule init
% git submodule update
```
There's more documentation on submodules here: https://git-scm.com/book/en/v2/Git-Tools-Submodules

Building the code requires you to have builds of the RDKit and Boost::Python available; the easiest way to get these is to set up a conda environment:

```
% conda create -n py312_shape python=3.12 cmake rdkit libboost-devel libboost-python-devel boost-cpp
% conda activate py312_shape
% wget -O $CONDA_PREFIX/include/rdkit/RDBoost/pyint_api.h https://raw.githubusercontent.com/rdkit/rdkit/Release_2024_03/Code/RDBoost/pyint_api.h
% mkdir build
% cd build
% cmake ..
% make -j4
% ctest --output-on-failure
```
**Note** that the horrible bit above with the `wget` should not be necessary after the 2024.03.6 release of the RDKit is available in conda-forge.

