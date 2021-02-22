# CalcLigRMSD: Calculate RMSD between prealigned molecules

`CalcLigRMSD` calculates the Root-mean-square deviation (RMSD) between two prealigned molecules. 

To evaluate the goodness of docking programs, docked protein-ligand structures are compared to the corresponding crystallographic complexes.
First, the protein structure is aligned based on the backbone atoms. Then, the RMSD between the docked and crystal ligands is calculated.
However, two problems may arise in the RMSD calculation, which are dealt with the `CalcLigRMSD` function: 

1. the atom names in the compared structures do not match, i.e. the same atoms present different names in the two coordinate files. In this case, we cannot rely on the atom names for the RMSD caculation and we need to first identify equivalent atoms. This "matching" problem is particularly challenging when the molecule of interest contains symmetrical chemical groups or when the entire molecule is symmetrical. In this case, the symmetry of the mocule needs to be taked into account. For example for a tri-methyl group, there are 6 possible symmetrical combinations: 
(C1, C2, C3), (C1, C3, C2), (C2, C1, C3), (C2, C3, C1), (C3, C1, C2), (C3, C2, C1). The `CalcLigRMSD` function takes symmetry into consideration (using rdkit functions) by calculating the RMSD with all possible "symmetrical" combinations and returns the minimum RMSD. You can see examples in the attached jupyter notebook Examples_CalcLigRMSD.ipynb.

2. one or both structures have missing atoms. This occurs, for example, in crystallographic structures when some atoms are not well defined in the electron density map. In this case, the `CalcLigRMSD` returns the RMSD for the maximum common substructure (MCS). You can see examples in the attached jupyter notebook Examples_CalcLigRMSD.ipynb.

### Dependencies
- rdkit
- numpy

### Usage
```
CalcLigRMSD(lig1, lig2, rename_lig2=True, output_filename='tmp.pdb')

Parameters
    ----------
    lig1 : RDKit molecule
    lig2 : RDKit molecule
    rename_lig2 : bool, optional
        True to rename the atoms of lig2 according to the atom names of lig1
    output_filename : str, optional
        If rename_lig2 is set to True, a PDB file with the renamed lig2 atoms is written as output.
        This may be useful to check that the RMSD has been "properly" calculated, 
        i.e. that the atoms have been properly matched for the calculation of the RMSD.
    
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two input molecules
```

### Example
```
from CalcLigRMSD import *

# load ligand pre-aligned coordinates
docked_ligand = Chem.rdmolfiles.MolFromPDBFile("data/docked_2c6e_JVE_pH74_netcharge1.pdb")
crystal_ligand = Chem.rdmolfiles.MolFromPDBFile("data/4uzh_JVE.pdb")
# calculate RMSD
ligs_rmsd = CalcLigRMSD(docked_ligand, crystal_ligand, rename_lig2 = True, output_filename = 'data/JVE_renamed.pdb')
print(f"RMSD: {ligs_rmsd:.2f}")

>> RMSD: 1.66
```

### Examples_CalcLigRMSD.ipynb
This jupyter notebook shows examples of how to use the `CalcLigRMSD` function.
In addition, it contains a python/pymol function that can be used to align protein structures and extract (pre-aligned) ligand coordinates.

### Acknowledgements
This work was conducted by Carmen Esposito as a member of the [Riniker lab](https://github.com/rinikerlab) at ETH Zurich.

### References
1. Velázquez-Libera, José Luis, et al. "LigRMSD: a web server for automatic structure matching and RMSD calculations among identical and similar compounds in protein-ligand docking." Bioinformatics 36.9 (2020): 2912-2914. [Link](https://academic.oup.com/bioinformatics/article/36/9/2912/5700716?login=true)
2. Bell, Eric W., and Yang Zhang. "DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism." Journal of cheminformatics 11.1 (2019): 1-9. [Link](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0362-7)
