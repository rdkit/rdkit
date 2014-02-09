An overview of the RDKit
%%%%%%%%%%%%%%%%%%%%%%%%

What is it?
===========

- Open source toolkit for cheminformatics

  - BSD licensed
  - Core data structures and algorithms in C++
  - Python (2.x) wrapper generated using Boost.Python
  - Java and C# wrappers generated with SWIG
  - 2D and 3D molecular operations
  - Descriptor generation for machine learning
  - Molecular database cartridge for PostgreSQL
  - Cheminformatics nodes for KNIME (distributed from the KNIME community site: http://tech.knime.org/community/rdkit)

- Operational:

  - http://www.rdkit.org
  - Supports Mac/Windows/Linux
  - Releases every 6 months
  - Web presence:

    - Homepage: http://www.rdkit.org
      
      Documentation, links

    - Github (https://github.com/rdkit)
 
      Bug tracker, git repository

    - Sourceforge (http://sourceforge.net/projects/rdkit) 
      
      Mailing lists, Downloads

    - Google code (http://code.google.com/p/rdkit/)
      
      wiki

  - Mailing lists at https://sourceforge.net/p/rdkit/mailman/, searchable archives available for
      `rdkit-discuss <http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/>`_ and
      `rdkit-devel <http://www.mail-archive.com/rdkit-devel@lists.sourceforge.net/>`_
       

- History:

  - 2000-2006: Developed and used at Rational Discovery for building predictive models for ADME, Tox, biological activity
  - June 2006: Open-source (BSD license) release of software, Rational Discovery shuts down
  - to present: Open-source development continues, use within Novartis, contributions from Novartis back to open-source version

Functionality overview
======================

- Input/Output: SMILES/SMARTS, SDF, TDT, SLN [1]_, Corina mol2 [1]_, PDB
- “Cheminformatics”:

  - Substructure searching
  - Canonical SMILES
  - Chirality support (i.e. R/S or E/Z labeling)
  - Chemical transformations (e.g. remove matching substructures)
  - Chemical reactions
  - Molecular serialization (e.g. mol <-> text)

- 2D depiction, including constrained depiction
- 2D->3D conversion/conformational analysis via distance geometry
- UFF and MMFF94/MMFF94S implementations for cleaning up structures
- Fingerprinting: Daylight-like, atom pairs, topological torsions, Morgan algorithm, “MACCS keys”, etc.
- Similarity/diversity picking
- 2D pharmacophores [1]_
- Gasteiger-Marsili charges
- Hierarchical subgraph/fragment analysis
- Bemis and Murcko scaffold determination
- RECAP and BRICS implementations
- Multi-molecule maximum common substructure [2]_
- Feature maps
- Shape-based similarity
- RMSD-based molecule-molecule alignment
- Shape-based alignment (subshape alignment [3]_) [1]_
- Unsupervised molecule-molecule alignment using Open3DAlign algorithm [4]_
- Integration with PyMOL for 3D visualization
- Functional group filtering
- Salt stripping
- Molecular descriptor library:

  - Topological (κ3, Balaban J, etc.)
  - Compositional (Number of Rings, Number of Aromatic Heterocycles, etc.)
  - Electrotopological state (Estate)
  - clogP, MR (Wildman and Crippen approach)
  - “MOE like” VSA descriptors
  - Feature-map vectors [5]_
  - MQN [6]_
- Similarity Maps [7]_

- Machine Learning:

  - Clustering (hierarchical, Butina)
  - Information theory (Shannon entropy, information gain, etc.)

- Tight integration with the `IPython <http://ipython.org>`_ notebook and qtconsole.


.. [1] These implementations are functional but are not necessarily the best, fastest, or most complete.

.. [2] Contribution from Andrew Dalke

.. [3] Putta, S., Eksterowicz, J., Lemmen, C. & Stanton, R. "A Novel Subshape Molecular Descriptor" *Journal of Chemical Information and Computer Sciences* **43:1623–35** (2003).

.. [4] Tosco, P., Balle, T. & Shiri, F. Open3DALIGN: an open-source software aimed at unsupervised ligand alignment. *J Comput Aided Mol Des* **25:777–83** (2011).

.. [5] Landrum, G., Penzotti, J. & Putta, S. "Feature-map vectors: a new class of informative descriptors for computational drug discovery" *Journal of Computer-Aided Molecular Design* **20:751–62** (2006).

.. [6] Nguyen, K. T., Blum, L. C., van Deursen, R. & Reymond, J.-L. Classification of Organic Molecules by Molecular Quantum Numbers. *ChemMedChem* **4:1803–5** (2009).

.. [7] Riniker, S. & Landrum, G. A. Similarity maps - a visualization strategy for molecular fingerprints and machine-learning methods. *Journal of Cheminformatics* **5:43** (2013).

The Contrib Directory
=====================

The Contrib directory, part of the standard RDKit distribution, includes code that has been contributed by members of the community.

- **LEF**: Local Environment Fingerprints 

  Contains python source code from the publications:

  - A. Vulpetti, U. Hommel, G. Landrum, R. Lewis and C. Dalvit, "Design and NMR-based screening of LEF, a library of chemical fragments with different Local Environment of Fluorine" *J. Am. Chem. Soc.* **131** (2009) 12949-12959. http://dx.doi.org/10.1021/ja905207t
  - A. Vulpetti, G. Landrum, S. Ruedisser, P. Erbel and C. Dalvit, "19F NMR Chemical Shift Prediction with Fluorine Fingerprint Descriptor" *J. of Fluorine Chemistry* **131** (2010) 570-577. http://dx.doi.org/10.1016/j.jfluchem.2009.12.024

  Contribution from Anna Vulpetti
  
- **M_Kossner**:

  Contains a set of pharmacophoric feature definitions as well as code for finding molecular frameworks.

  Contribution from Markus Kossner

- **PBF**: Plane of best fit

  Contains C++ source code and sample data from the publication: 

  N. C. Firth, N. Brown, and J. Blagg, "Plane of Best Fit: A Novel Method to Characterize the Three-Dimensionality of Molecules" *Journal of Chemical Information and Modeling* **52** 2516-2525 (2012). http://pubs.acs.org/doi/abs/10.1021/ci300293f

  Contribution from Nicholas Firth

- **mmpa**: Matched molecular pairs

  Python source and sample data for an implementation of the matched-molecular pair algorithm described in the publication:

  Hussain, J., & Rea, C. "Computationally efficient algorithm to identify matched molecular pairs (MMPs) in large data sets." *Journal of chemical information and modeling* **50** 339-348 (2010). http://dx.doi.org/10.1021/ci900450m

  Includes a fragment indexing algorithm from the publication:

  Wagener, M., & Lommerse, J. P. "The quest for bioisosteric replacements." *Journal of chemical information and modeling* **46** 677-685 (2006).

  Contribution from Jameed Hussain. 

- **SA_Score**: Synthetic assessibility score

  Python source for an implementation of the SA score algorithm described in the publication:

  Ertl, P. and Schuffenhauer A. "Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions" *Journal of Cheminformatics* **1:8** (2009)

  Contribution from Peter Ertl

- **fraggle**: A fragment-based molecular similarity algorithm

  Python source for an implementation of the fraggle similarity algorithm developed at GSK and described in this RDKit UGM presentation:
  https://github.com/rdkit/UGM_2013/blob/master/Presentations/Hussain.Fraggle.pdf

  Contribution from Jameed Hussain

  

License
=======

This document is copyright (C) 2013-2014 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself.
In simple words: “Do whatever you want with it, but please give us some credit.”
