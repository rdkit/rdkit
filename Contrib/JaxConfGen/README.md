GPU Conformer Generator
-----------------------

This is a highly experimental GPU-powered conformer generator, powered by jax/jaxlib, originally written by Yutong Zhao originally as part of the timemachine package. Functional forms may or may not be periodically updated. 

See requirements.txt on installation dependencies, and remember to initialize the git submodule:

git submodule init

To enable GPU usage, see:

https://github.com/google/jax#pip-installation

Though the JIT emitted CPU code is quite fast too.

Usage
-----

See demo.py on usage

Technical details
-----------------

The method proceeds by the classic Jeff Blaney DG method, where a symmetric NxN matrix is diagonalized, whose D dominant eigenvalues and eigenvectors are used to reconstruct an NxD starting conformer. 

This is then subsequently minimized by the SMIRNOFF Forcefield via the Adam optimizer.

Caveats
-------

1. There is a fairly expensive "up front" cost to emit the JIT kernels, this can be on the order of 1-2 seconds but is relatively insensitive to the size of the molecules or the batch size. So this is designed for large molecules and *lots* of conformers.

2. This is dependent on the SMIRNOFF99 forcefield developed by the Open Forcefield Consortium. While they do periodically update the toolkit (and parameters) we currently pin the submodule to a rather old version of the OFF toolkit (mainly because Yutong is too busy to deal with upgrades right now). 