# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Jameed Hussain, May 2013

import sys
from optparse import OptionParser
from rdkit import Chem
from rdkit import DataStructs
from collections import defaultdict

#input format
#query_substructs,query_smiles,SMILES,ID,Tversky_sim

#algorithm
#read in query_substructs and smiles
#feed to atomcontrib function to return generalised_SMILES
#use Tanimoto to compare generalised_SMILES with query smiles to give fraggle similarity

#atomcontrib algorithm
#generate fp of query_substructs (qfp)
#
#loop through atoms of smiles
#   For each atom
#	Generate partial fp of the atom (pfp)
#	Find Tversky sim of pfp in qfp
#	If Tversky < 0.8, mark atom in smiles
#
#Loop thru marked atoms
#	If marked atom in ring - turn all atoms in that ring to * (aromatic) or Sc (aliphatic)
#	For each marked atom
#		If aromatic turn to a *
#		If aliphatic turn to a Sc 
#
# Return modified smiles

def atomContrib(subs,smi,options):

    marked = {}

    def partialFP(atomID):

        #create empty fp
        modifiedFP = DataStructs.ExplicitBitVect(1024)

        modifiedFP.SetBitsFromList(aBits[atomID])

        tverskySim = DataStructs.TverskySimilarity(subsFp,modifiedFP,0,1)
 
        if(tverskySim < options.pfp):
            #print "%i %s: %f" % (atomID+1, pMol.GetAtomWithIdx(atomID).GetSymbol(), tverskySim)
            marked[atomID] = 1

    #generate mol object & fp for input mol
    aBits = [];

    pMol = Chem.MolFromSmiles(smi)
    pMolFp = Chem.RDKFingerprint(pMol, maxPath=5, fpSize=1024, nBitsPerHash=2, atomBits=aBits)

    #generate fp of query_substructs
    qsMol = Chem.MolFromSmiles(subs)
    subsFp = Chem.RDKFingerprint(qsMol, maxPath=5, fpSize=1024, nBitsPerHash=2)

    #loop through atoms of smiles and mark
    for atom in pMol.GetAtoms():           
        #store atoms to change
        partialFP(atom.GetIdx())

    #get rings to change
    ssr = pMol.GetRingInfo().AtomRings()

    #loop thru rings and records rings to change
    ringsToChange = {}
    for ringList in range(len(ssr)):
        #print "New ring"
        for ringAtom in range(len(ssr[ringList])):
            #print ssr[ringList][ringAtoms]
            if(  marked.has_key(ssr[ringList][ringAtom])  ):
                #print ssr[ringList][ringAtoms]
                ringsToChange[ringList] = 1

    #now add these ring atoms to marked
    for ringList in ringsToChange:
        for ringAtom in range(len(ssr[ringList])):
            marked[ ssr[ringList][ringAtom] ] = 1

    if(len(marked) > 0):
        #now mutate the marked atoms
        for key in marked:
            #print key
            if( pMol.GetAtomWithIdx(key).GetIsAromatic() ):
                #pMol.GetAtomWithIdx(key).SetAtomicNum(91)
                #this works everytime and causes far fewer problems
                pMol.GetAtomWithIdx(key).SetAtomicNum(0)
                pMol.GetAtomWithIdx(key).SetNoImplicit(True)
            else:
                #gives best sim
                pMol.GetAtomWithIdx(key).SetAtomicNum(21)
                #works better but when replace S it fails due to valency
                #pMol.GetAtomWithIdx(key).SetAtomicNum(6)

        try:
            Chem.SanitizeMol(pMol)
        except:
            sys.stderr.write("Can't parse smiles: %s\n" % (Chem.MolToSmiles(pMol)))
            pMol = Chem.MolFromSmiles(smi)

    return Chem.MolToSmiles(pMol)
			

parser = OptionParser(description="Program to post-process Tversky search results as part of Fraggle",
                    epilog="Format of input file: query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim\t"
                    "Output: SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity")
parser.add_option('-c','--cutoff',action='store', dest='cutoff', type='float', default=0.7,
                  help="Cutoff for fraggle similarity. Only results with similarity greater than the cutoff will be output. DEFAULT = 0.7")
parser.add_option('-p','--pfp',action='store', dest='pfp', type='float', default=0.8,
                  help="Cutoff for partial fp similarity. DEFAULT = 0.8")

if __name__ == '__main__':
    #parse the command line options
    (options, args) = parser.parse_args()

    if( (options.cutoff >= 0) and (options.cutoff <= 1) ):
        fraggle_cutoff = options.cutoff
    else:
        print "Fraggle cutoff must be in range 0-1"
        sys.exit(1)


    print "SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity"

    #create some data structure to store results
    id_to_smi = {}
    modified_query_fp = {}
    day_sim = {}
    frag_sim = {}
    query_size = {}
    query_mols = {}

    #generate dummy mol object which generates empty fp
    emptyMol = Chem.MolFromSmiles('*')

    #read the STDIN
    for line in sys.stdin:
        line = line.rstrip()
        qSubs,qSmi,qID,inSmi,id_,tversky = line.split(",")

        #add query to id_to_smi
        id_to_smi[qID] = qSmi
        id_to_smi[id_] = inSmi

        #add query to data structures
        frag_sim.setdefault(qID, defaultdict(float))
        day_sim.setdefault(qID, {})

        if(qID not in query_size):
            qMol = Chem.MolFromSmiles(qSmi)
            if(qMol == None):
                sys.stderr.write("Can't generate mol for: %s\n" % (qSmi) )
                continue
            query_mols[qID] = qMol
            query_size[qID] = qMol.GetNumAtoms()

        iMol = Chem.MolFromSmiles(inSmi)

        if(iMol == None):
            sys.stderr.write("Can't generate mol for: %s\n" % (inSmi) )
            continue

        #discard based on atom size
        if(iMol.GetNumAtoms() < query_size[qID]-3):
            #sys.stderr.write("Too small: %s\n" % (inSmi) )
            continue;

        if(iMol.GetNumAtoms() > query_size[qID]+4):
            #sys.stderr.write("Too large: %s\n" % (inSmi) )
            continue;
        
        #rdkit_sim,fraggle_sim = compute_fraggle_similarity()
        qFP = Chem.RDKFingerprint(query_mols[qID], maxPath=5, fpSize=1024, nBitsPerHash=2)
        iFP = Chem.RDKFingerprint(iMol, maxPath=5, fpSize=1024, nBitsPerHash=2)

        rdkit_sim = DataStructs.TanimotoSimilarity(qFP,iFP)

        #print "%s %s %s %s %f" % (qSmi,query_id,inSmi,id,rdkit_sim)
        #add to day_sim
        day_sim[qID][id_] = rdkit_sim

        #check if you have the fp for the modified query
        #and generate if need to
        qm_key = "%s_%s" % (qSubs,qSmi)
        if qm_key in modified_query_fp:
            qmMolFp = modified_query_fp[qm_key]
        else:
            query_modified = atomContrib(qSubs,qSmi,options)
            qmMol = Chem.MolFromSmiles(query_modified)
            #qmMolFp = FingerprintMols.FingerprintMol(qmMol)
            qmMolFp = Chem.RDKFingerprint(qmMol, maxPath=5, fpSize=1024, nBitsPerHash=2)
            #add to modified_query_fp
            modified_query_fp[qm_key] = qmMolFp

        retrieved_modified = atomContrib(qSubs,inSmi,options)
        #print "%s.%s" % (query_modified,retrieved_modified)
        rmMol = Chem.MolFromSmiles(retrieved_modified)

        #wrap in a try, catch
        try:
            rmMolFp = Chem.RDKFingerprint(rmMol, maxPath=5, fpSize=1024, nBitsPerHash=2)
            #print "%s," % (retrieved_modified),

            fraggle_sim=max(DataStructs.FingerprintSimilarity(qmMolFp,rmMolFp),
                            rdkit_sim)
                            
            #store fraggle sim if its the highest
            frag_sim[qID][id_] = max(frag_sim[qID][id_],fraggle_sim)
        except:
            sys.stderr.write("Can't generate fp for: %s\n" % (retrieved_modified) )

    #right, print out the results for the query
    #Format: SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,Daylight_Similarity
    for qID in frag_sim:
        for id_ in frag_sim[qID]:
            if(frag_sim[qID][id_] >= fraggle_cutoff):
                print "%s,%s,%s,%s,%s,%s" % (id_to_smi[id_],id_,id_to_smi[qID],qID,frag_sim[qID][id_],day_sim[qID][id_])

