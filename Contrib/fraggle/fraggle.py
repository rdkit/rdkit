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
import re
from rdkit import Chem

if (len(sys.argv) >= 2):
    print "Program to run the first part of Fraggle. Program splits the molecule\nready for the search\n";
    print "USAGE: ./fraggle.py <file_of_smiles";
    print "Format of smiles file: SMILES ID (space or comma separated)";
    print "Output: whole mol smiles,ID,fraggle split smiles\n";
    sys.exit(1)

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
def delete_bonds(bonds,type):

    #use the same parent mol object and create editable mol
    em = Chem.EditableMol(mol)

    #loop through the bonds to delete
    #print "Breaking bonds between atoms: ",bonds

    for i in bonds:
        #remove the bond
	em.RemoveBond(i[0],i[1])

	#now add attachement points
	newAtomA = em.AddAtom(Chem.Atom(0))
	em.AddBond(i[0],newAtomA,Chem.BondType.SINGLE)

	newAtomB = em.AddAtom(Chem.Atom(0))
	em.AddBond(i[1],newAtomB,Chem.BondType.SINGLE)

    #should be able to get away without sanitising mol
    #as the valencies should be okay
    modifiedMol = em.GetMol()

    #do not sanitise!
    #Chem.SanitizeMol(modifiedMol)

    fragmented_smi = Chem.MolToSmiles(modifiedMol,True)

    #print fragmented_smi
    fraggle_framentation = select_fragments(fragmented_smi,type)

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
    

def select_fragments(f_smi,type):

    result = ""
    result_hcount = 0
    fragments = f_smi.split(".")

    if(type == "acyclic"):
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

    elif(type == "cyclic"):
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


    elif(type == "cyclic_and_acyclic"):
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

#read the STDIN
for line in sys.stdin:

    line = line.rstrip()
    smi,id = re.split('\s|,',line)
    #print smi,id

    mol = Chem.MolFromSmiles(smi)

    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % (smi) )
	continue

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
    #print matching_atoms

    total_acyclic = len(acyclic_matching_atoms)
    bonds_selected = []

    #loop to generate every single and double cut in the molecule
    for x in xrange( total_acyclic ):
        #single cuts are not required
        #relevant single cut fragments can be found from the double cuts
        #for explanation see check_fragments method
        for y in xrange(x+1,total_acyclic):
            #print matching_atoms[x],matching_atoms[y]
            bonds_selected.append(acyclic_matching_atoms[x])
            bonds_selected.append(acyclic_matching_atoms[y])
            fragment = delete_bonds(bonds_selected,"acyclic")
            if fragment is not None:
                #print "%s" % (fragment)
                out_fragments.add(fragment)
            bonds_selected = []

    ##################################
    # Fused/Spiro exocyclic bond Cuts
    ##################################

    #find the relevant bonds to break
    cyclic_matching_atoms = mol.GetSubstructMatches(cyc_smarts)
    #print "Matching Atoms:"
    #print matching_atoms

    total_cyclic = len(cyclic_matching_atoms)
    bonds_selected = []

    #loop to generate every double cut of relevant bonds
    for x in xrange( total_cyclic ):
        for y in xrange(x+1,total_cyclic):
            #print matching_atoms[x],matching_atoms[y]
            bonds_selected.append(cyclic_matching_atoms[x])
            bonds_selected.append(cyclic_matching_atoms[y])
            fragment = delete_bonds(bonds_selected,"cyclic")
            bonds_selected = []

            if fragment is not None:
                #print "%s" % (fragment)
                out_fragments.add(fragment)

                #now do an acyclic cut with the successful cyclic cut
                for z in xrange(total_acyclic):
                    bonds_selected.append(cyclic_matching_atoms[x])
                    bonds_selected.append(cyclic_matching_atoms[y])
                    bonds_selected.append(acyclic_matching_atoms[z])
                    fragment = delete_bonds(bonds_selected,"cyclic_and_acyclic")
                    if fragment is not None:
                        #print "%s" % (fragment)
                        out_fragments.add(fragment)
                    bonds_selected = []

    #print out the unique fragments
    for x in out_fragments:
        #cansmi
        temp = Chem.MolFromSmiles(x)

        print "%s,%s,%s" % (smi,id,Chem.MolToSmiles(temp))

