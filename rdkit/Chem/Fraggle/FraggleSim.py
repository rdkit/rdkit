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
from rdkit import Chem,DataStructs

# our default rdkit fingerprinter parameters:
rdkitFpParams=dict(maxPath=5,fpSize=1024,nBitsPerHash=2)


#Fragmentation algorithm
#-----------------------

#identify acyclic bonds
#enumerate all single cuts
    #make sure you chop off more that 1 atom
    #keeps bits which are >60% query mol
#enumerate all double cuts
    #keeps bits with 1 attachment point (i.e throw middle bit away)
    #need to be >60% query mol

#identify exocyclic bonds
#enumerate all single "ring" cuts
    #Check if it results in more that one component
    #keep correct bit if >40% query mol

#enumerate successful "rings" cuts with an acyclic cut
    #Check if it results in more that one component
    #keep correct if >60% query mol

#start
def delete_bonds(mol,bonds,ftype,hac):

    #use the same parent mol object and create editable mol
    em = Chem.EditableMol(mol)

    #loop through the bonds to delete
    #print "Breaking bonds between atoms: ",bonds

    for b in bonds:
        #remove the bond
        em.RemoveBond(b[0],b[1])

        #now add attachement points
        newAtomA = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[0],newAtomA,Chem.BondType.SINGLE)

        newAtomB = em.AddAtom(Chem.Atom(0))
        em.AddBond(b[1],newAtomB,Chem.BondType.SINGLE)

    #should be able to get away without sanitising mol
    #as the valencies should be okay
    modifiedMol = em.GetMol()

    #do not do a full sanitization, but do find rings and calculate valences:
    Chem.SanitizeMol(modifiedMol,Chem.SanitizeFlags.SANITIZE_PROPERTIES|Chem.SanitizeFlags.SANITIZE_SYMMRINGS)

    fragmented_smi = Chem.MolToSmiles(modifiedMol,True)

    #print fragmented_smi
    fraggle_framentation = select_fragments(fragmented_smi,ftype,hac)

    return fraggle_framentation

def is_ring_cut_valid(smi):
    #to check is a fragment is a valid ring cut, it needs to match the
    #smarts: [$([#0][r].[r][#0]),$([#0][r][#0])]
    atom_count = 0
    valid = False

    m = Chem.MolFromSmiles(smi)
    if m is not None:
        #use gobal smarts
        if(m.HasSubstructMatch(cSma1) or m.HasSubstructMatch(cSma2)):
            atom_count = m.GetNumAtoms()
            valid = True

    return valid,atom_count
    

def select_fragments(f_smi,ftype,hac):

    result = ""
    result_hcount = 0
    fragments = f_smi.split(".")

    if(ftype == "acyclic"):
        for f in fragments:
            attachments = f.count("*")

            #check if terminal fragment
            if(attachments == 1):
                fMol = Chem.MolFromSmiles(f)
                fhac = fMol.GetNumAtoms()

                #if the fragment is 2 atoms (or less - includes attachement) it is too small
                #to be interesting. This check has the additional benefit
                #of pulling out the relevant single cuts as it discards
                #fragments where we only chop off a small part of the input cmpd
                if(fhac > 3):
                    result = "%s.%s" % (result,f)
                    result_hcount = result_hcount + fhac

        #needs to be greater than 60% of parent mol
        if( (result != "") and (result_hcount > 0.6*hac) ):
            #remove first character as it will always be "."
            result = result[1:]
        else:
            result = None

    elif(ftype == "cyclic"):
        result = None
        #make sure it is 2 components
        if( len(fragments) == 2):
            for f in fragments:
                #check if a valid cut
                valid,result_hcount = is_ring_cut_valid(f)

                if(valid):
                    #needs to be greater 3 heavy atoms and greater than 40% of parent mol
                    if((result_hcount > 3) and (result_hcount > 0.4*hac)):
                        result = f


    elif(ftype == "cyclic_and_acyclic"):
        #print f_smi
        result = ""

        #need to find the fragments which are valid which means they must be:
        #  Terminal (one attachment point) or valid ring cut
        for f in fragments:
            attachments = f.count("*")

            if(attachments >= 3):
                continue

            if(attachments == 2):
                #check if a valid cut
                valid,result_hcount = is_ring_cut_valid(f)
                if(valid):
                    #needs to be greater 3 heavy atoms
                    if(result_hcount > 3):
                        result = "%s.%s" % (result,f)

            elif(attachments == 1):
                fMol = Chem.MolFromSmiles(f)
                fhac = fMol.GetNumAtoms()
                #needs to be greater 3 heavy atoms
                if(fhac > 3):
                    result = "%s.%s" % (result,f)
                    result_hcount = result_hcount + fhac

        #print "F: %s" % (result)

        #appropriate fragmentations must have 2 components
        #result will always start with . because of the way it is constructed
        #hence 2 component result wil contain 2 dots
        if( (result != "") and (result.count(".") == 2) ):
            #take off the starting dot when building smiles
            fMol = Chem.MolFromSmiles(result[1:])
            result_hcount = fMol.GetNumAtoms()
            #needs to be greater 3 heavy atoms and greater than 60% of parent mol
            if((result_hcount > 3) and (result_hcount > 0.6*hac)):
                #take off the starting dot
                result = result[1:]
            else:
                result = None
        else:
            result = None

    return result

#Global smarts used by the program
#acyclic bond smarts
acyc_smarts = Chem.MolFromSmarts("[*]!@!=!#[*]")

#exocyclic/fused exocyclic bond smarts
cyc_smarts = Chem.MolFromSmarts("[R1,R2]@[r;!R1]")

#smarts used to find appropriate fragment for
#would use SMARTS: [$([#0][r].[r][#0]),$([#0][r][#0])]
#but rdkit doesn't support component SMARTS in recursive one - $([#0][r].[r][#0])
#hence split into two
cSma1 = Chem.MolFromSmarts("[#0][r].[r][#0]")
cSma2 = Chem.MolFromSmarts("[#0][r][#0]")

def generate_fraggle_fragmentation(mol):
    """
    >>> q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    >>> list(generate_fraggle_fragmentation(q))
     ['[*]C(=O)NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1', '[*]C(=O)c1cncc(C)c1.[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1', '[*]C(=O)c1cncc(C)c1.[*]Cc1cc(OC)c2ccccc2c1OC', '[*]C(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC', '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1', '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1', '[*]Cc1cc(OC)c2ccccc2c1OC.[*]NC(=O)c1cncc(C)c1', '[*]Cc1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1', '[*]N1CCC(NC(=O)c2cncc(C)c2)CC1.[*]c1cc(OC)c2ccccc2c1OC', '[*]NC(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC', '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1', '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1', '[*]c1c(CN2CCC(NC(=O)c3cncc(C)c3)CC2)cc(OC)c2ccccc12', '[*]c1c(OC)cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c1[*]', '[*]c1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12', '[*]c1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1']
    """
    #query mol heavy atom count
    hac = mol.GetNumAtoms()

    #different cuts can give the same fragments
    #to use out_fragments to remove them
    out_fragments = set()

    ######################
    # Single acyclic Cuts
    ######################

    #find the relevant bonds to break
    acyclic_matching_atoms = mol.GetSubstructMatches(acyc_smarts)
    #print "Matching Atoms:"
    #print("acyclic matching atoms: ",acyclic_matching_atoms)

    total_acyclic = len(acyclic_matching_atoms)
    bonds_selected = []

    #loop to generate every single and double cut in the molecule
    for x in range( total_acyclic ):
        #single cuts are not required
        #relevant single cut fragments can be found from the double cuts
        #for explanation see check_fragments method
        for y in range(x+1,total_acyclic):
            #print matching_atoms[x],matching_atoms[y]
            bonds_selected.append(acyclic_matching_atoms[x])
            bonds_selected.append(acyclic_matching_atoms[y])
            fragment = delete_bonds(mol,bonds_selected,"acyclic",hac)
            if fragment is not None:
                #print(fragment)
                out_fragments.add(fragment)
            bonds_selected = []

    #print(out_fragments)
    ##################################
    # Fused/Spiro exocyclic bond Cuts
    ##################################

    #find the relevant bonds to break
    cyclic_matching_atoms = mol.GetSubstructMatches(cyc_smarts)
    #print("cyclic matching atoms: ",cyclic_matching_atoms)
    #print "Matching Atoms:"
    #print matching_atoms

    total_cyclic = len(cyclic_matching_atoms)
    bonds_selected = []

    #loop to generate every double cut of relevant bonds
    for x in range( total_cyclic ):
        for y in range(x+1,total_cyclic):
            #print matching_atoms[x],matching_atoms[y]
            bonds_selected.append(cyclic_matching_atoms[x])
            bonds_selected.append(cyclic_matching_atoms[y])
            fragment = delete_bonds(mol,bonds_selected,"cyclic",hac)
            bonds_selected = []

            if fragment is not None:
                #print "%s" % (fragment)
                out_fragments.add(fragment)

                #now do an acyclic cut with the successful cyclic cut
                for z in range(total_acyclic):
                    bonds_selected.append(cyclic_matching_atoms[x])
                    bonds_selected.append(cyclic_matching_atoms[y])
                    bonds_selected.append(acyclic_matching_atoms[z])
                    fragment = delete_bonds(mol,bonds_selected,"cyclic_and_acyclic",hac)
                    if fragment is not None:
                        #print "%s" % (fragment)
                        out_fragments.add(fragment)
                    bonds_selected = []

    return sorted(list(out_fragments))

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
def atomContrib(subs,mol,tverskyThresh=0.8):
    marked = {}

    def partialFP(atomID,tverskyThresh):

        #create empty fp
        modifiedFP = DataStructs.ExplicitBitVect(1024)

        modifiedFP.SetBitsFromList(aBits[atomID])

        tverskySim = DataStructs.TverskySimilarity(subsFp,modifiedFP,0,1)
 
        if(tverskySim < tverskyThresh):
            #print "%i %s: %f" % (atomID+1, pMol.GetAtomWithIdx(atomID).GetSymbol(), tverskySim)
            marked[atomID] = 1

    #generate mol object & fp for input mol
    aBits = [];

    pMol = Chem.Mol(mol.ToBinary())
    pMolFp = Chem.RDKFingerprint(pMol, atomBits=aBits, **rdkitFpParams)

    #generate fp of query_substructs
    qsMol = Chem.MolFromSmiles(subs)
    subsFp = Chem.RDKFingerprint(qsMol, **rdkitFpParams)

    #loop through atoms of smiles and mark
    for atom in pMol.GetAtoms():           
        #store atoms to change
        partialFP(atom.GetIdx(),tverskyThresh)

    #get rings to change
    ssr = pMol.GetRingInfo().AtomRings()

    #loop thru rings and records rings to change
    ringsToChange = {}
    for ringList in range(len(ssr)):
        #print "New ring"
        for ringAtom in range(len(ssr[ringList])):
            #print ssr[ringList][ringAtoms]
            if ssr[ringList][ringAtom] in marked:
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
        except Exception:
            sys.stderr.write("Can't parse smiles: %s\n" % (Chem.MolToSmiles(pMol)))
            pMol = Chem.Mol(mol.ToBinary())

    return pMol


modified_query_fps = {}
def compute_fraggle_similarity_for_subs(inMol,qMol,qSmi,qSubs,tverskyThresh=0.8):
    qFP = Chem.RDKFingerprint(qMol, **rdkitFpParams)
    iFP = Chem.RDKFingerprint(inMol, **rdkitFpParams)

    rdkit_sim = DataStructs.TanimotoSimilarity(qFP,iFP)

    qm_key = "%s_%s" % (qSubs,qSmi)
    if qm_key in modified_query_fps:
        qmMolFp = modified_query_fps[qm_key]
    else:
        qmMol  = atomContrib(qSubs,qMol,tverskyThresh)
        qmMolFp = Chem.RDKFingerprint(qmMol, **rdkitFpParams)
        modified_query_fps[qm_key] = qmMolFp

    rmMol = atomContrib(qSubs,inMol,tverskyThresh)

    #wrap in a try, catch
    try:
        rmMolFp = Chem.RDKFingerprint(rmMol, **rdkitFpParams)

        fraggle_sim=max(DataStructs.FingerprintSimilarity(qmMolFp,rmMolFp),
                        rdkit_sim)
        #print '\t',qSubs,fraggle_sim,rdkit_sim

    except Exception:
        sys.stderr.write("Can't generate fp for: %s\n" % (Chem.MolToSmiles(rmMol,True)))
        fraggle_sim = 0.0

    return rdkit_sim,fraggle_sim

def GetFraggleSimilarity(queryMol,refMol,tverskyThresh=0.8):
    """ return the Fraggle similarity between two molecules

    >>> q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    >>> m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.980...
    >>> match
    '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1'    

    >>> m = Chem.MolFromSmiles('COc1cc(CN2CCC(Nc3nc4ccccc4s3)CC2)c(OC)c2ccccc12')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.794...
    >>> match
    '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1'    

    >>> q = Chem.MolFromSmiles('COc1ccccc1')
    >>> sim,match = GetFraggleSimilarity(q,m)
    >>> sim
    0.347...
    >>> match
    '[*]c1ccccc1'    
 
    """
    if hasattr(queryMol,'_fraggleDecomp'):
        frags = queryMol._fraggleDecomp
    else:
        frags = generate_fraggle_fragmentation(queryMol)
        queryMol._fraggleDecomp = frags
    qSmi = Chem.MolToSmiles(queryMol,True)
    result=0.0
    bestMatch=None
    for frag in frags:
        rdksim,fragsim= compute_fraggle_similarity_for_subs(refMol,queryMol,qSmi,frag,tverskyThresh)
        if fragsim>result:
            result=fragsim
            bestMatch=frag
    return result,bestMatch
        
    

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"],optionflags=doctest.ELLIPSIS+doctest.NORMALIZE_WHITESPACE)


if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)

