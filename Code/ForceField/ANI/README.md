ANI Style Neural Networks as a ForceField
================================
[Manan Goel](manan.goel@research.iiit.ac.in)

This work was done as a part of Google Summer of Code 2020.
<br>

Introduction
--------------
[ANI](https://pubs.rsc.org/en/content/articlelanding/2017/sc/c6sc05720a#!divAbstract) uses Neural Network potentials found using a highly modified version of the Behler and Parrinello symmetry functions to build single-atom Atomic Environment Vectors. This work adds support for the [ANI-1x and ANI-1ccx](https://chemrxiv.org/articles/The_ANI-1ccx_and_ANI-1x_Data_Sets_Coupled-Cluster_and_Density_Functional_Theory_Properties_for_Molecules/10050737) architectures which are ensembles of 8 models trained on atoms with H, C, N and O only. It has been shown through that outputs from ANI are chemically accurate compared in reference to DFT calculations while also providing a major speed up.
<br>

Implementation Details
---------------
This work is an extension the ForceField API of RDKit which had implementations of MMFF and UFF. As the aforementioned Force Fields have different kinds of contributions, for ANI the total energy is split into individual atomic contributions which are then summed together to get the total energy.
<br>
The parameters for the ANI-1x and ANI-1ccx neural networks as well as the AEV generation are obtained from the open source implementation of ANI in pytorch - [torchani](https://github.com/aiqm/torchani)
<br>
All the linear algebra is done using the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library including storing the Atomic Environment Vectors and the Neural Network weights.
<br>
This involves the following steps - 
- Generating the Atomic Environment Vectors from the atomic coordinates and atom types in the molecule according to the parameters for the Gaussian functions(Code can be found in ```Code/GraphMol/Descriptors/AtomicEnvironmentVector.cpp``` and ```Code/GraphMol/Descriptors/AtomicEnvironmentVector.h```).
- Storing and Loading the weights for the Neural Network for which the Boost Serialization Library is used (Code can be found in ```Code/Numerics/EigenSerializer```).
- Forward Propogation of the obtained AEVs through the neural network which is simple matrix multiplication(Code can be found in ```Code/ForceField/ANI``` and ```Code/GraphMol/ForceFieldHelpers/ANI```).
- Finding derivatives of the obtained energy in order to find the force on each atom to get the lowest energy conformers(Code can be found in ```Code/ForceField/ANI``` and ```Code/GraphMol/ForceFieldHelpers/ANI```).
- Implemeting python wrappers for each of the functionalities.

The obtained energy and gradient calculation functions had to be implemented since the ForceField API already was connected to the BFGS optimizer.
<br>

Advantages
----------------
- Adds a new accurate way of finding energy of molecules as well as optimizing the structure of atoms with H, C, N and O.
- Adds an interface for integrating further more accurate forms of ANI style neural networks like the recently published ANI-2x.
<br>

Scopes of Improvement
----------------
- The AEV generation code can still be heavily optimized to make it faster. This section of the pipeline is a significant bottlneck.
- As of now, the derivative calculation is numeric i.e. each derivative calculation requires 6 energy calculations which is extremely time consuming. This can be replaced with using analytical derivatives which would make it faster.
- The ```EigenSerializer``` module is not supported on windows and hence ANI cannot be built on windows.

