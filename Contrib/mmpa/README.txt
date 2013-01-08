This directory contains the scripts used to generate matched molecular pairs (MMPs) from an input list of SMILES. 
The fragment indexing algorithm used in the scripts is described in the following publications:

Hussain, J., & Rea, C. (2010). "Computationally efficient algorithm to identify matched molecular pairs (MMPs) 
in large data sets." Journal of chemical information and modeling, 50(3), 339-348.
http://dx.doi.org/10.1021/ci900450m
	
Wagener, M., & Lommerse, J. P. (2006). "The quest for bioisosteric replacements." 
Journal of chemical information and modeling, 46(2), 677-685.

The scripts requires RDKit (www.rdkit.org) be installed and properly configured.

Help is available for all the scripts using the -h option

To find all the MMPs in your set 
--------------------------------

The program to generate the MMPs from a set is divided into two parts; fragmentation and indexing. 

Before running the programs, make sure your input set of SMILES:
    - does not contain mixtures (salts etc.) 
    - does not contain "*" atoms
    - has been canonicalised using RDKit.

If your smiles set doesn't satisfy the conditions above the programs are likely to fail or in the case of 
canonicalisation result in not identifying MMPs involving H atom substitution.

1) Fragmentation command:

python rfrag.py <SMILES_FILE >FRAGMENT_OUTPUT

Example command: 
python rfrag.py <sample.smi >sample_fragmented.txt

Format of SMILES_FILE is: SMILES ID <space or comma separated>
See sample.smi for an example input file

Format of output: WHOLE_MOL_SMILES,ID,SMILES_OF_CORE,SMILES_OF_CONTEXT
See sample_fragmented.txt for an example output file 

2) Index command:

python indexing.py <FRAGMENT_OUTPUT >MMP_OUTPUT.CSV

Format of output:
SMILES_OF_LEFT_MMP,SMILES_OF_RIGHT_MMP,ID_OF_LEFT_MMP,ID_OF_RIGHT_MMP,SMIRKS_OF_TRANSFORMATION,SMILES_OF_CONTEXT 

This program has several options (see help from program below):

Usage: indexing.py [options]

Program to generate MMPs

Options:
  -h, --help            show this help message and exit
  -s, --symmetric       Output symmetrically equivalent MMPs, i.e output both
                        cmpd1,cmpd2, SMIRKS:A>>B and cmpd2,cmpd1, SMIRKS:B>>A
  -m MAXSIZE, --maxsize=MAXSIZE
                        Maximum size of change (in heavy atoms) allowed in
                        matched molecular pairs identified. DEFAULT=10.
                        Note: This option overrides the ratio option if both
                        are specified.
  -r RATIO, --ratio=RATIO
                        Maximum ratio of change allowed in matched molecular
                        pairs identified. The ratio is: size of change /
                        size of cmpd (in terms of heavy atoms). DEFAULT=0.3.
                        Note: If this option is used with the maxsize option,
                        the maxsize option will be used.

Example commands (with sample outputs):

Default settings:
python indexing.py <sample_fragmented.txt >sample_mmps_default.csv

Output symmetrically equivalent MMPs (ie forward and reverse transforms):
python indexing.py -s <sample_fragmented.txt >sample_mmps_sym.csv

Output MMPs where maximum size of change is 3 heavy atoms:
python indexing.py -m 3 <sample_fragmented.txt >sample_mmps_maxheavy.csv

Output MMPs where no more that 10% of the compound has changed:
python indexing.py -r 0.1 <sample_fragmented.txt >sample_mmps_maxratio.csv

Output symmetrically equivalent MMPs where maximum size of change is 3 heavy atoms:
python indexing.py -s -m 3 <sample_fragmented.txt >sample_mmps_sym_maxheavy.csv
 
SMIRKS canonicalisation
-----------------------

The MMP identification script uses a SMIRKS canonicalisation routine so the same change always has the same output SMIRKS.

To canonicalise a SMIRKS (generated elsewhere) so it is in the same format as MMP identification scripts use command:

python cansmirks.py <SMIRKS_FILE >SMIRKS_OUTPUT_FILE

Example command:
python cansmirk.py <sample_smirks.txt >sample_cansmirks.txt

Format of SMIRKS_FILE: SMIRKS ID <space or comma separated>
See sample_smirks.txt for an example input file
  
Format of output: CANONICALISED_SMIRKS ID
See sample_cansmirks.txt for an example output file

Note: The script will NOT deal with SMARTS characters, so the SMIRKS must contain valid SMILES for left and right hand sides.

The algorithm used to canonicalise SMIRKS is as follows:
1) Canonicalise the LHS.
2) For the LHS the 1st asterisk (attachment point) in the SMILES will have label 1, 2nd asterisk will have label 2 and so on
3) For the RHS, if you have a choice (ie. two attachement points are symmetrically equivalent), always put the label 
with lower numerical value on the earlier attachement point in the canonicalised SMILES

  
In the event you use the scripts for publication please reference the original publication:

Hussain, J., & Rea, C. (2010). "Computationally efficient algorithm to identify matched molecular pairs (MMPs) 
in large data sets." Journal of chemical information and modeling, 50(3), 339-348.
http://dx.doi.org/10.1021/ci900450m
