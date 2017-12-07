import argparse
import subprocess
import pickle
import os
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser( "Fast clustering for chemoinformatics" )
parser.add_argument( "input" )
parser.add_argument( "n", help = "the number of clusters" )
#parser.add_argument( "-l", help = "limit value of cluster bisection" )
#parser.add_argument( "-p", help = "output similarity points" )
parser.add_argument( "--output", default = "cluster.tsv" )
parser.add_argument( "-c", help = "filename of centroid", default = "centroid.tsv" )

def smi2fp( molid, smiles ):
    mol = Chem.MolFromSmiles( smiles )
    onbits = AllChem.GetMorganFingerprintAsBitVect( mol, 2 ).GetOnBits()
    row = molid
    for bit in onbits:
        row += "\tFP_{}\t1.0".format( bit )
    row += "\n"
    return row


if __name__ == "__main__":
    currentdir = os.getcwd()
    args = parser.parse_args()
    inf = open( args.input, "r" )
    n = args.n
    c = args.c
    output = args.output
    tempf = open( "fp.tsv", "w" )
    for line in inf:
        line = line.rstrip().split( "\t" )
        tempf.write( smi2fp( line[0], line[1] ))
    tempf.close()
    res = subprocess.call( "time bayon -p -c {} -n {} fp.tsv > {}".format( c, n, output ), shell=True )

    #parse results
    parsefile = open( output.split(".")[0]+"_parse.tsv", "w" )
    inputf = open( output, "r" )
    for line in inputf:
        line = line.rstrip().split( "\t" )
        cluster_id = line[0]
        for i in range( 1, len( line )-1, 2 ) :
            molid = line[ i ]
            point = line[ i + 1 ]
            #print(itr)
            parsefile.write("{}\t{}\tCLS_{}\n".format( molid, point, cluster_id ))
    parsefile.close()


    if res == 0:
        print( "Done!" )
    else:
        print( "Error :-(" )
