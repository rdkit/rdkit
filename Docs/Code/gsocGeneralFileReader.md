**Do Not Review**

# [GSoC 2020] General File Reader and Multithreaded Mol Supplier
## Overview
The General File Reader, as the name suggests, provides the user with the appropriate `MolSupplier` object to parse a file of a given format. Thus for instance earlier if one wanted to parse a file of smiles, say `input.smi`, then one would need to explicitly construct an object `SmilesMolSupplier`. However, with the implementation provided in the GeneralFileReader, one can easily pass the file name along with supplier options to obtain the appropriate `MolSupplier` object determined by the file format. Furthermore, the General File Reader also provides an interface with the `MultithreadedMolSupplier` objects for Smiles and SDF file formats. Besides the implementation, test cases are also included to demonstrate the utility of the General File Reader.

The Multithreaded Mol Supplier provides a concurrent implementation of the usual base class `MolSupplier`. Due to time constraints, multithreaded versions of only Smiles, and SD Mol Suppliers were implemented. The motivation for this part stemmed from parsing large Smiles or SDF files. With the current implementation the user, for instance, can construct the object `MultithreadedSmilesMolSupplier` to parse a smiles file with a large number of records. Besides the implementation, test cases are also included to demonstrate the correctness and performance of the `MultithreadedMolSupplier`. Here is a brief summary of the performance result obtained by running the function `testPerformance` on @greglandrum's machine:

```
Duration for SmilesMolSupplier: 6256 (milliseconds)
Maximum Duration for MultithreadedSmilesMolSupplier: 6972 (milliseconds) with 1 writer thread
Minimum Duration for MultithreadedSmilesMolSupplier: 855 (milliseconds) with 15 writer threads

Duration for SDMolSupplier: 2584 (milliseconds) 
Maximum Duration for MultithreadedSDMolSupplier: 2784 (milliseconds) with 1 writer thread
Minimum Duration for MultithreadedSDMolSupplier: 729 (milliseconds) with 7 writer threads
```

## Implementation
Implementation of the General File Reader is quite concise and makes use of only two methods `determineFormat` and `getSupplier`. The former determines the file and the compression format given pathname, while the latter returns a pointer to `MolSupplier` object given pathname and `SupplierOptions`.

Regarding the implementation of the `MultithreadedMolSupplier`, the first step was to implement a thread-safe blocking queue of fixed capacity. This would allow us to extract and process records from the file concurrently. The concurrent queue was implemented with a single lock and two condition variables to signal whether the queue was empty or full. Test cases checking the correctness of the `ConcurrentQueue` are also included in the project.

The next step required the implementation of the base class `MultithreadedMolSupplier` which would manage the input and output queue. The input queue would be populated by the method `extractNextRecord` that would read a record from a given file/stream, whereas the output queue would be populated by the method `processMoleculeRecord` that would first pop a record from the input queue and then process it into an object of type `ROMol *`. The reader thread would thus call `extractNextRecord` until no record can be read, while the writer thread(s) would call the method `processMoleculeRecord` until the output queue is done and empty. The child classes `MultithreadedSmilesMolSupplier` and `MultithreadedSDMolSupplier` primarily provide implementations of the methods, `extractNextRecord` and `processMoleculeRecord`. Both suppliers were tested on various files with different parameter values for input queue size, output queue size, and the number of writer threads. 

## Further Work
Due to time constraints and the difficulty involved in debugging concurrent code, there were a few things that could not be completed.
1. In cases where the file format is less defined, it might be useful to parse the file content to discover the file format and possible Supplier options. The current implementation does not support this and only uses the pathname to determine the appropriate Supplier.
2. Wrappers for the Multithreaded Smiles and SD Suppliers in other languages such as Java were not implemented in this project.

## Changes made for the General File Reader and Multithreaded Mol Supplier:

List of important files added:

- [GeneralFileReader.h](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/GeneralFileReader.h) and [testGeneralFileReader.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/testGeneralFileReader.cpp)
- [ConcurrentQueue.h](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/RDGeneral/ConcurrentQueue.h) and [testConcurrentQueue.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/RDGeneral/testConcurrentQueue.cpp)
- [MultithreadedMolSupplier.h](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedMolSupplier.h) and [MultithreadedMolSupplier.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedMolSupplier.cpp)
- [MultithreadedSmilesMolSupplier.h](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedSmilesMolSupplier.h), [MultithreadedSmilesMolSupplier.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedSmilesMolSupplier.cpp) and its [Python wrapper](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/Wrap/MultithreadedSmilesMolSupplier.cpp)
- [MultithreadedSDMolSupplier.h](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.h), [MultithreadedSDMolSupplier.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp) and its [Python wrapper](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/Wrap/MultithreadedSDMolSupplier.cpp).
- [testMultithreadedMolSupplier.cpp](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/FileParsers/testMultithreadedMolSupplier.cpp) and [testMultithreadedMolSupplier.py](https://github.com/shrey183/rdkit/blob/GSOC-2020/Code/GraphMol/Wrap/testMultithreadedMolSupplier.py) for testing Python wrappers.



#### Any other comments?
This project is still work in progress.

