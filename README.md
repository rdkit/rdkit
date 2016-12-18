# RDKit
[![Build status](https://travis-ci.org/rdkit/rdkit.svg)](https://travis-ci.org/rdkit/rdkit)
[![Documentation Status](https://readthedocs.org/projects/rdkit/badge/?version=latest)](http://rdkit.readthedocs.org/en/latest/)

[RDKit](https://github.com/rdkit/rdkit) is a collection of cheminformatics and machine-learning software written in C++ and Python.

  * [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt) - a business friendly license for open source
  * Core data structures and algorithms in C++
  * [Python (2.x and 3.x) wrapper](http://www.rdkit.org/docs/GettingStartedInPython.html) generated using Boost.Python
  * Java and C# wrappers generated with SWIG
  * 2D and 3D molecular operations
  * [Descriptor](http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors) and [Fingerprint](http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints) generation for machine learning
  * Molecular database [cartridge](http://www.rdkit.org/docs/Cartridge.html) for PostgreSQL supporting substructure and similarity searches as well as many descriptor calculators
  * Cheminformatics nodes for [KNIME](http://tech.knime.org/community/rdkit)
  * [Contrib](https://github.com/rdkit/rdkit/tree/master/Contrib) folder with useful community-contributed software harnessing the power of the RDKit

## Web presence

  * [RDKit web page](https://github.com/rdkit/rdkit)
  * [Blog](https://rdkit.blogspot.com)
  * [Documentation](http://www.rdkit.org/docs/index.html) also
    available at [ReadTheDocs](http://rdkit.readthedocs.org/en/latest/)

#### Code

  * [GitHub code](https://github.com/rdkit) and [bug tracker](https://github.com/rdkit/rdkit/issues)

#### Community

  * [Mailing list](https://sourceforge.net/p/rdkit/mailman/)
  * [Google +](https://plus.google.com/u/0/116996224395614252219)
  * [LinkedIn](https://www.linkedin.com/groups/RDKit-8192558/about)

#### Materials from user group meetings

  * [2012 UGM](http://www.rdkit.org/UGM/2012/)
  * [2013 UGM](https://github.com/rdkit/UGM_2013)
  * [2014 UGM](https://github.com/rdkit/UGM_2014)
  * [2015 UGM](https://github.com/rdkit/UGM_2015)
  * [2016 UGM](https://github.com/rdkit/UGM_2016)

## Documentation
Available on the [RDKit page](http://www.rdkit.org/docs/index.html)
and in the [Docs](https://github.com/rdkit/rdkit/tree/master/Docs) folder on GitHub

## Installation

Installation instructions are available in [Docs/Book/Install.md](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md).

#### Binary distributions, anaconda, homebrew

  * Windows binaries are available with each [release](https://github.com/rdkit/rdkit/releases).
  * [RPMs](https://copr.fedoraproject.org/coprs/giallu/rdkit/) for RedHat Enterprise Linux, Centos, and Fedora. Contributed by Gianluca Sforna.
  * [homebrew](https://github.com/rdkit/homebrew-rdkit) formula for building on the Mac. Contributed by Eddie Cao.
  * [recipes](https://github.com/rdkit/conda-rdkit) for building using the excellent conda package manager. Contributed by Riccardo Vianello.

## Projects using RDKit

  * [ChEMBL Beaker](https://github.com/mnowotka/chembl_beaker) - standalone web server wrapper for RDKit and OSRA
  * [myChEMBL](https://github.com/chembl/mychembl) ([blog post](http://chembl.blogspot.de/2013/10/chembl-virtual-machine-aka-mychembl.html), [paper](http://bioinformatics.oxfordjournals.org/content/early/2013/11/20/bioinformatics.btt666)) - A virtual machine implementation of open data and cheminformatics tools
  * [sdf_viewer.py](https://github.com/apahl/sdf_viewer) - an interactive SDF viewer
  * [sdf2ppt](https://github.com/dkuhn/sdf2ppt) - Reads an SDFile and displays molecules as image grid in powerpoint/openoffice presentation.
  * [MolGears](https://github.com/admed/molgears) - A cheminformatics tool for bioactive molecules
  * [PYPL](http://www.biochemfusion.com/downloads/#OracleUtilities) - Simple cartridge that lets you call Python scripts from Oracle PL/SQL.
  * [shape-it-rdkit](https://github.com/jandom/shape-it-rdkit) - Gaussian molecular overlap code shape-it (from silicos it) ported to RDKit backend
  * [WONKA](http://wonka.sgc.ox.ac.uk/WONKA/) - Tool for analysis and interrogation of protein-ligand crystal structures
  * [OOMMPPAA](http://oommppaa.sgc.ox.ac.uk/OOMMPPAA/) - Tool for directed synthesis and data analysis based on protein-ligand crystal structures


## License

Code released under the [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt).
