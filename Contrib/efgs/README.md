Source code of EFGs, a complete an accurate implementation of Ertl's functional group (FG) detection algorithm.

For a RDKit molecule, it provides four outputs: 
a) a PNG binary string with an image of the molecule with color-highlighted
functional groups; 
b) a list of sets of atom indices (idx), each set corresponding to a
functional group; 
c) a list of pseudo-SMILES canonicalized strings for the full functional
groups; 
d) a list of RDKit labeled mol objects, one for each full functional group 

efgs.py: code of get_dec_fgs function and additional auxiliary functions

try_efgs.py: example code to execute get_dec_fgs in Python with a Pandas dataframe

ch33query.sql: sql query used to extract ChEMBL33 bioactive compounds

ChemRxiv preprint for full details on the implementation: 
EFGs: A Complete and Accurate Implementation of Ertl’s Functional Group Detection Algorithm in RDKit
https://chemrxiv.org/engage/chemrxiv/article-details/67543afe7be152b1d03a820c
