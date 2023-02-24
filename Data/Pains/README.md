# Original data source:
Saubern, S., Guha, R. & Baell, J. B. KNIME Workflow to Assess PAINS
Filters in SMARTS Format. Comparison of RDKit and Indigo
Cheminformatics Libraries. Mol. Inform. 30, 847–850 (2011). 

  The paper: http://onlinelibrary.wiley.com/doi/10.1002/minf.201100076/full
  The workflow: http://www.myexperiment.org/workflows/1841.html

# Modifications
Many of the individual SMARTS definitions have been edited in order to
get them working properly with the RDKit.

Note that these SMARTS contain explicit hydrogen atoms, so are best
used with merged H queries, such as with the mergeHs=True option in
MolFromSmarts.  See Greg's blogpost
http://rdkit.blogspot.com/2015/08/curating-pains-filters.html for more
information.
