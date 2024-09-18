# RDKit
[![Azure build Status](https://dev.azure.com/rdkit-builds/RDKit/_apis/build/status/rdkit.rdkit?branchName=master)](https://dev.azure.com/rdkit-builds/RDKit/_build/latest?definitionId=1&branchName=master)
[![DOI](https://zenodo.org/badge/10009991.svg)](https://zenodo.org/badge/latestdoi/10009991)


## What is it?

The [RDKit](https://www.rdkit.org) is a collection of cheminformatics and machine-learning software written in C++ and Python.

  * [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt) - a business friendly license for open source
  * Core data structures and algorithms in C++
  * [Python 3.x wrapper](https://www.rdkit.org/docs/GettingStartedInPython.html) generated using Boost.Python
  * Java and C# wrappers generated with SWIG
  * JavaScript and CFFI wrappers around important functionality
  * 2D and 3D molecular operations
  * [Descriptor](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors) and [Fingerprint](http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints) generation for machine learning
  * Molecular database [cartridge](https://www.rdkit.org/docs/Cartridge.html) for PostgreSQL supporting substructure and similarity searches as well as many descriptor calculators
  * Cheminformatics nodes for [KNIME](https://www.knime.com/rdkit)
  * [Contrib](https://github.com/rdkit/rdkit/tree/master/Contrib) folder with useful community-contributed software harnessing the power of the RDKit


## Installation and getting started

If you are working in Python and using conda (our recommendation), installation is super easy:

```shell-session
$ conda install -c conda-forge rdkit
```

You can then take a look at our [Getting Started in Python](https://rdkit.org/docs/GettingStartedInPython.html) guide.

More detailed installation instructions are available in [Docs/Book/Install.md](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md).

## Documentation
Available on the [RDKit page](https://www.rdkit.org/docs/index.html)
and in the [Docs](https://github.com/rdkit/rdkit/tree/master/Docs) folder on GitHub

The [RDKit blog](https://greglandrum.github.io/rdkit-blog/) often has useful tips and tricks.

## Support and Community

If you have questions, comments, or suggestions, the best places for those are:

  * The [mailing list](https://sourceforge.net/p/rdkit/mailman/)
  * [GitHub discussions](https://github.com/rdkit/rdkit/discussions)

If you've found a bug or would like to request a feature, please [create an issue](https://github.com/rdkit/rdkit/issues)

We also have a [LinkedIn group](https://www.linkedin.com/groups/RDKit-8192558/about)


## License

Code released under the [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt).
