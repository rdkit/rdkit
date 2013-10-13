Fraggle Readme
--------------

This directory contains the scripts used to run Fraggle.

The algorithm used in the scripts was described at the 2nd RDKit UGM (October 
2013). The presentation can be found at:
https://github.com/rdkit/UGM_2013/blob/master/Presentations/Hussain.Fraggle.pdf

The benchmarking carried out in the presentation utilised the open source
benchmarking platform described in:

Riniker, Sereina, and Gregory A. Landrum. "Open-source platform to benchmark
fingerprints for ligand-based virtual screening." Journal of cheminformatics 5.1
(2013): 26.

With the addition of the following scripts:

fraggle.py
cxn_tversky.py
atomcontrib.py
calculate_scored_lists_mod.py

The information below describes how to run the Fraggle similarity algorithm with
a query compound against a file of database compounds.

How to run Fraggle:
-------------------

Fraggle works in three steps:

1) Need to fragment your query molecule(s)
2) Run a Tversky Search using the generated fragments
3) Post-process results of the Tversky search to give final output

It is recommended to run a standard RDK5 similarity alongside Fraggle

The scripts requires RDKit (www.rdkit.org) be installed and properly configured.
Help is available for all the scripts using the -h option

Step 1
------

Command:
python fraggle.py <QUERY_FILE >FRAGGLE_FRAGMENTS

Exmaple command:
python fraggle.py < data/query.smi > data/query_fragmentation.csv

Format of QUERY_FILE is: SMILES ID <space or comma separated>
See query.smi for an example input file

Format of FRAGGLE_FRAGMENTS: whole mol smiles,ID,fraggle split smiles
See query_fragmentation.csv for an example output file

The following help is available using the -h option:

Program to run the first part of Fraggle. Program splits the molecule
ready for the search

USAGE: ./fraggle.py <file_of_smiles
Format of smiles file: SMILES ID (space or comma separated)
Output: whole mol smiles,ID,fraggle split smiles

Step 2
------

The second step to take the fragments generated in step 1 and run a Tversky
search against your database of molecules. This is the rate determining step of
the algorithm so it is recommended to do this against a database with an
appropriate chemistry cartridge. However a python script is provided which
utilises RDKit.

The script uses a default tversky cut-off of 0.8 (alpha=0,beta=1) which seems to
work the reasonably well for the rdk5 fp.

Command:
python rdkit_tversky.py -f FRAGGLE_FRAGMENTS <DB_SMILES_FILE >TVERSKY_OUTPUT

Example command:
python rdkit_tversky.py -f data/query_fragmentation.csv < data/ChEMBL_11265_actives.smi > data/fragmentation_tversky_out

Format of FRAGGLE_FRAGMENTS file is: whole mol smiles,ID,fraggle split smiles
See query_fragmentation.csv for an example file

Format of DB_SMILES_FILE: SMILES ID (space or comma separated)
See ChEMBL_11265_actives.smi for an example file

Format of TVERSKY_OUTPUT: query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim
See fragmentation_tversky_out for an example file

The following help is available using the -h option:
Usage: rdkit_tversky.py [options]

Program to Tversky search results as part of Fraggle

Options:
  -h, --help            show this help message and exit
  -f F_FILE, --frags=F_FILE
                        File containing the query fragmentations from Fraggle
  -c CUTOFF, --cutoff=CUTOFF
                        Cutoff for Tversy similarity. Only Tversky results
                        with similarity greater than the cutoff will be
                        output. DEFAULT = 0.8

Format of input file: whole mol smiles,ID,fraggle split smiles
Output:
query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim

Step 3
------

The last step is to perform the post-processing to give you the final Fraggle
similarity

Command:
python atomcontrib.py <TVERSKY_OUTPUT >FINAL_FRAGGLE_RESULTS

Example command:
python atomcontrib.py < data/fragmentation_tversky_out > data/final_fraggle_results.csv

Format of TVERSKY_OUTPUT file is:
query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim
See fragmentation_tversky_out for an example file

Format of FINAL_FRAGGLE_RESULTS:
SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity
See final_fraggle_results.csv for an example file

This program has several options (see help from program below):

Usage: atomcontrib.py [options]

Program to post-process Tversky search results as part of Fraggle

Options:
  -h, --help            show this help message and exit
  -c CUTOFF, --cutoff=CUTOFF
                        Cutoff for fraggle similarity. Only results with
                        similarity greater than the cutoff will be output.
                        DEFAULT = 0.7
  -p PFP, --pfp=PFP     Cutoff for partial fp similarity. DEFAULT = 0.8

Format of input file:
query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim
Output: SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity

