# RDKit
[![Azure build Status](https://dev.azure.com/rdkit-builds/RDKit/_apis/build/status/rdkit.rdkit?branchName=master)](https://dev.azure.com/rdkit-builds/RDKit/_build/latest?definitionId=1&branchName=master)
[![Documentation Status](https://readthedocs.org/projects/rdkit/badge/?version=latest)](http://rdkit.readthedocs.org/en/latest/)
[![DOI](https://zenodo.org/badge/10009991.svg)](https://zenodo.org/badge/latestdoi/10009991)


[RDKit](https://github.com/rdkit/rdkit) is a collection of cheminformatics and machine-learning software written in C++ and Python.

  * [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt) - a business friendly license for open source
  * Core data structures and algorithms in C++
  * [Python 3.x wrapper](https://www.rdkit.org/docs/GettingStartedInPython.html) generated using Boost.Python
  * Java and C# wrappers generated with SWIG
  * 2D and 3D molecular operations
  * [Descriptor](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors) and [Fingerprint](http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints) generation for machine learning
  * Molecular database [cartridge](https://www.rdkit.org/docs/Cartridge.html) for PostgreSQL supporting substructure and similarity searches as well as many descriptor calculators
  * Cheminformatics nodes for [KNIME](https://www.knime.com/rdkit)
  * [Contrib](https://github.com/rdkit/rdkit/tree/master/Contrib) folder with useful community-contributed software harnessing the power of the RDKit

#### Community

  * [Mailing list](https://sourceforge.net/p/rdkit/mailman/)
  * [GitHub discussions](https://github.com/rdkit/rdkit/discussions)
  * [LinkedIn](https://www.linkedin.com/groups/RDKit-8192558/about)

#### Code

  * [GitHub code](https://github.com/rdkit) and [bug tracker](https://github.com/rdkit/rdkit/issues)

## Web presence

  * [RDKit web page](https://rdkit.org)
  * [Blog](https://rdkit.blogspot.com)
  * [Documentation](https://www.rdkit.org/docs/index.html)

#### Materials from user group meetings

  * [2012 UGM](http://www.rdkit.org/UGM/2012/)
  * [2013 UGM](https://github.com/rdkit/UGM_2013)
  * [2014 UGM](https://github.com/rdkit/UGM_2014)
  * [2015 UGM](https://github.com/rdkit/UGM_2015)
  * [2016 UGM](https://github.com/rdkit/UGM_2016)
  * [2017 UGM](https://github.com/rdkit/UGM_2017)
  * [2018 UGM](https://github.com/rdkit/UGM_2018)
  * [2019 UGM](https://github.com/rdkit/UGM_2019)
  * [2020 UGM](https://github.com/rdkit/UGM_2020)
  * [2021 UGM](https://github.com/rdkit/UGM_2021)

## Documentation
Available on the [RDKit page](https://www.rdkit.org/docs/index.html)
and in the [Docs](https://github.com/rdkit/rdkit/tree/master/Docs) folder on GitHub

## Installation

Installation instructions are available in [Docs/Book/Install.md](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md).

#### Binary distributions, anaconda, homebrew

  * [binaries for conda python](https://anaconda.org/rdkit/rdkit) or, if you are using the conda-forge stack, the RDKit is also [available from conda-forge](https://anaconda.org/conda-forge/rdkit).
  * [RPMs](https://src.fedoraproject.org/rpms/rdkit) for RedHat Enterprise Linux, Centos, and Fedora. Contributed by Gianluca Sforna.
  * [debs](https://blends.debian.org/debichem/tasks/cheminformatics) for Ubuntu and other Debian-derived Linux distros. Contributed by the Debichem team.
  * [homebrew](https://github.com/rdkit/homebrew-rdkit) formula for building on the Mac. Contributed by Eddie Cao.
  * [recipes](https://github.com/rdkit/conda-rdkit) for building using the excellent conda package manager. Contributed by Riccardo Vianello.
  * [APKs](https://github.com/daverona/alpine-rdkit) for Alpine Linux. Contributed by da Verona
  * [Wheels](https://pypi.org/project/rdkit/) at PyPi for all major platforms and python versions. Contributed by Christopher Kuenneth

## Projects using RDKit

- [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) - ChEMBL protocols used to standardise and salt strip molecules.
- [FPSim2](https://github.com/chembl/FPSim2) - Simple package for fast molecular similarity searches.
- [Datamol](https://datamol.io/) ([docs](https://doc.datamol.io/stable/), [repo](https://github.com/datamol-org/datamol/)) - A Python library to intuitively manipulate molecules.
- [Scopy](https://github.com/kotori-y/Scopy) ([docs](https://scopy.iamkotori.com/), [paper](https://doi.org/10.1093/bib/bbaa194)) - an integrated negative design Python library for desirable HTS/VS database design
- [stk](https://github.com/lukasturcani/stk) ([docs](https://stk.readthedocs.io), [paper](https://onlinelibrary.wiley.com/doi/10.1002/jcc.25377)) -
a Python library for building, manipulating, analyzing and automatic design of molecules.
- [gpusimilarity](https://github.com/schrodinger/gpusimilarity) - A Cuda/Thrust implementation of fingerprint similarity searching
- [Samson Connect](https://www.samson-connect.net) - Software for adaptive modeling and simulation of nanosystems
- [mol_frame](https://github.com/apahl/mol_frame) - Chemical Structure Handling for Dask and Pandas DataFrames
- [RDKit.js](https://github.com/rdkit/rdkit-js) - The official JavaScript release of RDKit
- [DeepChem](https://deepchem.io) - python library for deep learning for chemistry
- [mmpdb](https://github.com/rdkit/mmpdb) - Matched molecular pair database generation and analysis
- [CheTo](https://github.com/rdkit/CheTo) ([paper](http://pubs.acs.org/doi/10.1021/acs.jcim.7b00249))- Chemical topic modeling
- [OCEAN](https://github.com/rdkit/OCEAN) ([paper](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00067))- Optimized cross reactivity estimation
- [ChEMBL Beaker](https://github.com/mnowotka/chembl_beaker) - standalone web server wrapper for RDKit and OSRA
- [ZINC](http://zinc15.docking.org) - Free database of commercially-available compounds for virtual screening
- [sdf_viewer.py](https://github.com/apahl/sdf_viewer) - an interactive SDF viewer
- [sdf2ppt](https://github.com/dkuhn/sdf2ppt) - Reads an SDFile and displays molecules as image grid in powerpoint/openoffice presentation.
- [MolGears](https://github.com/admed/molgears) - A cheminformatics tool for bioactive molecules
- [PYPL](http://www.biochemfusion.com/downloads/#OracleUtilities) - Simple cartridge that lets you call Python scripts from Oracle PL/SQL.
- [shape-it-rdkit](https://github.com/jandom/shape-it-rdkit) - Gaussian molecular overlap code shape-it (from silicos it) ported to RDKit backend
- [WONKA](http://wonka.sgc.ox.ac.uk/WONKA/) - Tool for analysis and interrogation of protein-ligand crystal structures
- [OOMMPPAA](http://oommppaa.sgc.ox.ac.uk/OOMMPPAA/) - Tool for directed synthesis and data analysis based on protein-ligand crystal structures
- [OCEAN](https://github.com/rdkit/OCEAN) - web-tool for target-prediction of chemical structures which uses ChEMBL as datasource
- [chemfp](http://chemfp.com) - very fast fingerprint searching
- [rdkit_ipynb_tools](https://github.com/apahl/rdkit_ipynb_tools) - RDKit Tools for the IPython Notebook
- [Vernalis KNIME nodes](https://www.knime.com/book/vernalis-nodes-for-knime-trusted-extension)
- [Erlwood KNIME nodes](https://www.knime.com/community/erlwood)
- [AZOrange](https://github.com/AZcompTox/AZOrange)

## License

Code released under the [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt).
