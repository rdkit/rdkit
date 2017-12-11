/*
This script performs fast clustering of SMILES

Clustering method is repeated bi section, the method looks like -k-means.
To use this script, the user needs to install bayon at first.
Input format: Tab separated SMILES strings (SMILES \t molID \n ...)
Please see more details in README.

*/

import argparse
import subprocess
import pickle
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def getArgParser():
    """ Create the argument parser """
    parser = argparse.ArgumentParser("Fast clustering for chemoinformatics")
    parser.add_argument("input", help="filename of input file")
    parser.add_argument("nclusters", metavar="N", help="the number of clusters")
    parser.add_argument("--output",  help="file name of output, tab deliminated format", default="clustered.tsv")
    parser.add_argument("-c", "--centroid", metavar="CENTROID",  help="filename of centroid information. tab deliminated format", default="centroid.tsv")
    return parser

def smi2fp(molid, smiles):
    mol = Chem.MolFromSmiles(smiles)
    onbits = AllChem.GetMorganFingerprintAsBitVect(mol, 2).GetOnBits()
    row = molid
    for bit in onbits:
        row += "\tFP_{}\t1.0".format(bit)
    row += "\n"
    return row


if __name__ == "__main__":
    parser = getArgParser()
    args = parser.parse_args()
    inputf = open(args.input, "r")
    nclusters = args.nclusters
    centroid = args.centroid
    output = args.output
    tempf = open("fp.tsv", "w")
    for line in inputf:
        line = line.rstrip().split("\t")
        tempf.write(smi2fp(line[0], line[1]))
    tempf.close()
    res = subprocess.call("time bayon -p -c {} -n {} fp.tsv > {}".format(centroid, nclusters, output), shell=True)

    #parse results
    parsefile = open(output.split(".")[0]+"_parse.tsv", "w")
    inputf = open(output, "r")
    for line in inputf:
        line = line.rstrip().split("\t")
        cluster_id = line[0]
        for i in range(1, len(line)-1, 2) :
            molid = line[ i ]
            point = line[ i + 1 ]
            parsefile.write("{}\t{}\tCLS_ID_{}\n".format(molid, point, cluster_id))
    parsefile.close()


    if res == 0:
        print("Done!")
    else:
        print( "-----" )
