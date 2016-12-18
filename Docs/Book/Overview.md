# An overview of the RDKit

## What is it?

### Open source toolkit for cheminformatics
-   Business-friendly BSD license
-   Core data structures and algorithms in C++
-   Python (2.x and 3.x) wrapper generated using Boost.Python
-   Java and C\# wrappers generated with SWIG
-   2D and 3D molecular operations
-   Descriptor generation for machine learning
-   Molecular database cartridge for PostgreSQL
-   Cheminformatics nodes for KNIME (distributed from the KNIME community site: http://tech.knime.org/community/rdkit)

### Operational:
- http://www.rdkit.org
- Supports Mac/Windows/Linux
- Releases every 6 months
- Web presence:
    - Homepage: http://www.rdkit.org
      Documentation, links
    - Github (https://github.com/rdkit)
      Downloads, bug tracker, git repository
    - Sourceforge (http://sourceforge.net/projects/rdkit)
      Mailing lists
    - Blog (https://rdkit.blogspot.com)
      Tips, tricks, random stuff
    - Tutorials (https://github.com/rdkit/rdkit-tutorials)
      Jupyter-based tutorials for using the RDKit
    - KNIME integration (https://github.com/rdkit/knime-rdkit)
      RDKit nodes for KNIME
- Mailing lists at https://sourceforge.net/p/rdkit/mailman/, searchable archives available for [rdkit-discuss](http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/) and [rdkit-devel](http://www.mail-archive.com/rdkit-devel@lists.sourceforge.net/)
- Social media:
    - Twitter: @RDKit_org
    - LinkedIn: https://www.linkedin.com/groups/8192558
    - Google+: https://plus.google.com/u/0/116996224395614252219
    - Slack: https://rdkit.slack.com (invite required, contact Greg)

### History:
-   2000-2006: Developed and used at Rational Discovery for building predictive models for ADME, Tox, biological activity
-   June 2006: Open-source (BSD license) release of software, Rational Discovery shuts down
-   to present: Open-source development continues, use within Novartis, contributions from Novartis back to open-source version

## Functionality overview
## Basics
- Input/Output: SMILES/SMARTS, SDF, TDT, SLN [1](#footnote1), Corina mol2 [1](#footnote1), PDB, sequence notation, FASTA (peptides only), HELM (peptides only)
- Substructure searching
- Canonical SMILES
- Chirality support (i.e. R/S or E/Z labeling)
- Chemical transformations (e.g. remove matching substructures)
- Chemical reactions
- Molecular serialization (e.g. mol \<-\> text)
- 2D depiction, including constrained depiction
- Fingerprinting: Daylight-like, atom pairs, topological torsions, Morgan algorithm, “MACCS keys”, extended reduced graphs, etc.
- Similarity/diversity picking
- Gasteiger-Marsili charges
- Bemis and Murcko scaffold determination
- Salt stripping
- Functional-group filters

### 2D
- 2D pharmacophores [1](#footnote1)
- Hierarchical subgraph/fragment analysis
- RECAP and BRICS implementations
- Multi-molecule maximum common substructure [2](#footnote2)
- Enumeration of molecular resonance structures
- Molecular descriptor library:
  - Topological (κ3, Balaban J, etc.)
  - Compositional (Number of Rings, Number of Aromatic Heterocycles, etc.)
  - Electrotopological state (Estate)
  - clogP, MR (Wildman and Crippen approach)
  - “MOE like” VSA descriptors
  - MQN [6](#footnote6)
- Similarity Maps [7](#footnote7)
- Machine Learning:
  - Clustering (hierarchical, Butina)
  - Information theory (Shannon entropy, information gain, etc.)
- Tight integration with the [Jupyter](http://jupyter.org) notebook (formerly the IPython notebook) and [Pandas](http://pandas.pydata.org/).

### 3D
- 2D-\>3D conversion/conformational analysis via distance geometry, including optional use of experimental torsion angle potentials [9](#footnote9)
- UFF and MMFF94/MMFF94S implementations for cleaning up structures
- Pharmacophore embedding (generate a pose of a molecule that matches a 3D pharmacophore) [1](#footnote1)
- Feature maps
- Shape-based similarity
- RMSD-based molecule-molecule alignment
- Shape-based alignment (subshape alignment [3](#footnote3)) [1](#footnote1)
- Unsupervised molecule-molecule alignment using the Open3DAlign algorithm [4](#footnote4)
- Integration with PyMOL for 3D visualization
- Molecular descriptor library:
  - Moments-of-inertia based descriptors: PMI, NPR, PBF, etc.
  - Feature-map vectors [5](#footnote5)
- Torsion Fingerprint Differences for comparing conformations [8](#footnote8)

### Integration with other open-source projects
- [KNIME](https://tech.knime.org/community/rdkit): Workflow and analytics tool
- [Django](http://django-rdkit.readthedocs.org/en/latest/): "The web framework for perfectionists with deadlines"
- [PostgreSQL](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Cartridge.rst): Extensible relational database
- [Lucene](https://github.com/rdkit/org.rdkit.lucene): Text-search engine [1](#footnote1)

### Usage by other open-source projects
- [ChEMBL Beaker](https://github.com/mnowotka/chembl_beaker) - standalone web server wrapper for RDKit and OSRA
- [myChEMBL](https://github.com/chembl/mychembl) ([blog post](http://chembl.blogspot.de/2013/10/chembl-virtual-machine-aka-mychembl.html), [paper](http://bioinformatics.oxfordjournals.org/content/early/2013/11/20/bioinformatics.btt666)) - A virtual machine implementation of open data and cheminformatics tools
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
- [chemicalite](https://github.com/rvianello/chemicalite) - SQLite integration for the RDKit
- [Vernalis KNIME nodes](https://tech.knime.org/book/vernalis-nodes-for-knime-trusted-extension)
- [Erlwood KNIME nodes](https://tech.knime.org/community/erlwood)
- [AZOrange](https://github.com/AZcompTox/AZOrange)


## The Contrib Directory

The Contrib directory, part of the standard RDKit distribution, includes code that has been contributed by members of the community.

### LEF: Local Environment Fingerprints

Contains python source code from the publications:

-   A. Vulpetti, U. Hommel, G. Landrum, R. Lewis and C. Dalvit, "Design and NMR-based screening of LEF, a library of chemical fragments with different Local Environment of Fluorine" *J. Am. Chem. Soc.* **131** (2009) 12949-12959. http://dx.doi.org/10.1021/ja905207t
-   Vulpetti, G. Landrum, S. Ruedisser, P. Erbel and C. Dalvit, "19F NMR Chemical Shift Prediction with Fluorine Fingerprint Descriptor" *J. of Fluorine Chemistry* **131** (2010) 570-577. http://dx.doi.org/10.1016/j.jfluchem.2009.12.024

Contribution from Anna Vulpetti

### M\_Kossner

Contains a set of pharmacophoric feature definitions as well as code for finding molecular frameworks.

Contribution from Markus Kossner

### PBF: Plane of best fit

Contribution from Nicholas Firth

*Note* as of the 2016.09.1 release this functionality is part of the RDKit core.

Contains C++ source code and sample data from the publication:

Firth, N. Brown, and J. Blagg, "Plane of Best Fit: A Novel Method to Characterize the Three-Dimensionality of Molecules" *Journal of Chemical Information and Modeling* **52** 2516-2525 (2012). http://pubs.acs.org/doi/abs/10.1021/ci300293f


### mmpa: Matched molecular pairs

Python source and sample data for an implementation of the matched-molecular pair algorithm described in the publication:

Hussain, J., & Rea, C. "Computationally efficient algorithm to identify matched molecular pairs (MMPs) in large data sets." *Journal of chemical information and modeling* **50** 339-348 (2010). http://dx.doi.org/10.1021/ci900450m

Includes a fragment indexing algorithm from the publication:

Wagener, M., & Lommerse, J. P. "The quest for bioisosteric replacements." *Journal of chemical information and modeling* **46** 677-685 (2006).

Contribution from Jameed Hussain.

### SA\_Score: Synthetic assessibility score

Python source for an implementation of the SA score algorithm described in the publication:

Ertl, P. and Schuffenhauer A. "Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions" *Journal of Cheminformatics* **1:8** (2009)

Contribution from Peter Ertl

### fraggle: A fragment-based molecular similarity algorithm

Python source for an implementation of the fraggle similarity algorithm developed at GSK and described in this RDKit UGM presentation: https://github.com/rdkit/UGM_2013/blob/master/Presentations/Hussain.Fraggle.pdf

Contribution from Jameed Hussain

### pzc: Tools for building and validating classifiers

Contribution from Paul Czodrowski

### ConformerParser: parser for Amber trajectory files

Contribution from Sereina Riniker

*Note* as of the 2016.09.1 release this functionality is part of the RDKit core.

### NP_Score: Natural-product likeness score

Python source for an implementation of the NP score algorithm described in the publication:

"Natural Product Likeness Score and Its Application for Prioritization of Compound Libraries"
Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer
*Journal of Chemical Information and Modeling* **48:68-74** (2008)
http://pubs.acs.org/doi/abs/10.1021/ci700286x

Contribution from Peter Ertl

### AtomAtomSimilarity: atom-atom-path method for fragment similarity

Python source for an implementation of the Atom-Atom-Path similarity method for fragments described in the publication:

Gobbi, A., Giannetti, A. M., Chen, H. & Lee, M.-L. "Atom-Atom-Path similarity and Sphere Exclusion clustering: tools for prioritizing fragment hits." *J. Cheminformatics* **7** 11 (2015). http://dx.doi.org10.1186/s13321-015-0056-8

Contribution from Richard Hall

## Footnotes

<a name="footnote1">1</a>: These implementations are functional but are not necessarily the best, fastest, or most complete.

<a name="footnote2">2</a>: Originally contributed by Andrew Dalke

<a name="footnote3">3</a>: Putta, S., Eksterowicz, J., Lemmen, C. & Stanton, R. "A Novel Subshape Molecular Descriptor" *Journal of Chemical Information and Computer Sciences* **43:1623–35** (2003).

<a name="footnote4">4</a>: Tosco, P., Balle, T. & Shiri, F. "Open3DALIGN: an open-source software aimed at unsupervised ligand alignment." *J Comput Aided Mol Des* **25:777–83** (2011).

<a name="footnote5">5</a>: Landrum, G., Penzotti, J. & Putta, S. "Feature-map vectors: a new class of informative descriptors for computational drug discovery" *Journal of Computer-Aided Molecular Design* **20:751–62** (2006).

<a name="footnote6">6</a>: Nguyen, K. T., Blum, L. C., van Deursen, R. & Reymond, J.-L. "Classification of Organic Molecules by Molecular Quantum Numbers." *ChemMedChem* **4:1803–5** (2009).

<a name="footnote7">7</a>: Riniker, S. & Landrum, G. A. "Similarity maps - a visualization strategy for molecular fingerprints and machine-learning methods." *Journal of Cheminformatics* **5:43** (2013).

<a name="footnote8">8</a>: Schulz-Gasch, T., Schärfer, C., Guba, W. & Rarey, M. "TFD: Torsion Fingerprints As a New Measure To Compare Small Molecule Conformations." *J. Chem. Inf. Model.* **52:1499–1512** (2012).

<a name="footnote9">9</a>: Riniker, S. & Landrum, G. A. "Better informed distance geometry: Using what we know to improve conformation generation." *J. Chem. Inf. Model.* **55:2562–74** (2015).

## License

This document is copyright (C) 2013-2016 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/4.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
