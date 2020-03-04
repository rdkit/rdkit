# RDKit Coding Standards and Guidelines

## What is this?

This document is an overview of the suggested coding standards and conventions for the RDKit project. This isn't comprehensive and not all existing RDKit code follows all standards (ah, the joys of a long history).

If you find mistakes, or have suggestions for improvements, please either fix them yourselves in the source document (the .md file) or send them to the mailing list: <rdkit-devel@lists.sourceforge.net> (you will need to subscribe first)

## C++

### A note on Modern C++

In order to continue to allow the code to be built on older systems, the RDKit is currently using a subset of what's available in modern C++ compilers (by which I mean C++11, C++14, etc.). The expectation is that it should be possible to build the RDKit with g++ v4.8.

### Naming conventions

We use camelCase for everything.
Class and namespace names start with an upper-case letter, i.e. `ExplicitBitVect` or `MolFragmenter`.
Function and method names start with a lower-case letter, e.g. `ROMol::getNumAtoms()`.

### Formatting

We are following the formatting suggestions from [Google's C++ style guidelines](https://google.github.io/styleguide/cppguide.html). We strongly suggest using an auto-formatter like [clang-format](http://clang.llvm.org/docs/ClangFormat.html) to do the formatting. Clang-format can be integrated with editors like emacs, vim, or atom to automatically format code as you write/save it. It's an amazing time saver! *Note* that at least version 3.8 of clang-format is required, unfortunately the configuration files are not backwards compatible.

### Testing

### Use of language features

### Other specifics


## Python

### Naming conventions

### Formatting

### Testing

### Use of language features

### Other specifics


## License

This document is copyright (C) 2015 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/3.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
