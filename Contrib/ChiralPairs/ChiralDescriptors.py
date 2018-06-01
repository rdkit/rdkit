#
#  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
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
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
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
# Created by Nadine Schneider & Peter Ertl, July 2017


from collections import defaultdict, Counter, namedtuple

import seaborn as sns
import numpy as np
import re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# build an svg grid image to print
def _svgsToGrid(svgs, labels, svgsPerRow=4,molSize=(250,150),fontSize=12):
    
    matcher = re.compile(r'^(<.*>\n)(<rect .*</rect>\n)(.*)</svg>',re.DOTALL) 
    hdr='' 
    ftr='</svg>' 
    rect='' 
    nRows = len(svgs)//svgsPerRow 
    if len(svgs)%svgsPerRow : nRows+=1 
    blocks = ['']*(nRows*svgsPerRow)
    labelSizeDist = fontSize*5
    fullSize=(svgsPerRow*(molSize[0]+molSize[0]/10.0),nRows*(molSize[1]+labelSizeDist))
    print(fullSize)

    count=0
    for svg,name in zip(svgs,labels):
        h,r,b = matcher.match(svg).groups()
        if hdr == '': 
            hdr = h.replace("width='"+str(molSize[0])+"px'","width='%dpx'"%fullSize[0])
            hdr = hdr.replace("height='"+str(molSize[1])+"px'","height='%dpx'"%fullSize[1])
        if rect == '': 
            rect = r
        legend = '<text font-family="sans-serif" font-size="'+str(fontSize)+'px" text-anchor="middle" fill="black">\n'
        legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*2)+'">'+name.split('|')[0]+'</tspan>\n'
        if len(name.split('|')) > 1:
            legend += '<tspan x="'+str(molSize[0]/2.)+'" y="'+str(molSize[1]+fontSize*3.5)+'">'+name.split('|')[1]+'</tspan>\n'
        legend += '</text>\n'
        blocks[count] = b + legend
        count+=1

    for i,elem in enumerate(blocks): 
        row = i//svgsPerRow 
        col = i%svgsPerRow 
        elem = rect+elem 
        blocks[i] = '<g transform="translate(%d,%d)" >%s</g>'%(col*(molSize[0]+molSize[0]/10.0),row*(molSize[1]+labelSizeDist),elem) 
    res = hdr + '\n'.join(blocks)+ftr 
    return res 


def determineAtomSubstituents(atomID, mol, distanceMatrix, verbose=False):
    atomPaths = distanceMatrix[atomID]
    # determine the direct neighbors of the atom
    neighbors = [n for n,i in enumerate(atomPaths) if i == 1]
    # store the ids of the neighbors (substituents)
    sub_dict = defaultdict(list)
    # track in how many substituents an atom is involved (can happen in rings)
    sharedNeighbors_dict = defaultdict(int)
    # determine the max path lenght for each substituent
    maxShell_dict=defaultdict(int)
    for n in neighbors:
        sub_dict[n].append(n)
        sharedNeighbors_dict[n]+=1
        maxShell_dict[n]=0
    # second shell of neighbors
    mindist=2
    # max distance from atom
    maxdist=int(np.max(atomPaths))
    for d in range(mindist,maxdist+1):
        if verbose:
            print("Shell: ",d)
        newShell = [n for n,i in enumerate(atomPaths) if i == d]
        for aidx in newShell:
            if verbose:
                print("Atom ", aidx," in shell ",d)
            atom = mol.GetAtomWithIdx(aidx)
            # find neighbors of the current atom that are part of the substituent already
            for n in atom.GetNeighbors():
                nidx = n.GetIdx()
                for k,v in sub_dict.items():
                    # is the neighbor in the substituent and is not inthe same shell as the current atom 
                    # and we haven't added the current atom already then put it in the correct substituent list
                    if nidx in v and nidx not in newShell and aidx not in v:
                        sub_dict[k].append(aidx)
                        sharedNeighbors_dict[aidx]+=1
                        maxShell_dict[k]=d
                        if verbose:
                            print("Atom ",aidx," assigned to ",nidx)
    if verbose:
        print(sub_dict)
        print(sharedNeighbors_dict)
        
    return sub_dict, sharedNeighbors_dict, maxShell_dict

def _getSizeOfSubstituents(sub, sharedNeighbors_dict, weighdownShared=True):
    size=0
    if weighdownShared:
        for a in sub:
            size+=1.0/sharedNeighbors_dict[a]
        return size
    return(len(sub))
    
def getBondsSubstituent(mol, atoms):
    bonds=[]
    for b in mol.GetBonds():
        a1 = b.GetBeginAtomIdx()
        a2 = b.GetEndAtomIdx()
        if a1 in atoms and a2 in atoms:
            bonds.append(b.GetIdx())
    return bonds

def getAromaticBondsSubstituent(mol, subAtoms):
    bonds = getBondsSubstituent(mol, subAtoms)
    aroBonds=0
    for b in bonds:
        if mol.GetBondWithIdx(b).GetIsAromatic():
            aroBonds+=1
    return aroBonds

def getRotatableBondsSubstituent(mol, subAtoms):
    rotatableBond = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    matches = mol.GetSubstructMatches(rotatableBond)
    numRotBonds=0
    for a1,a2 in matches:
        if a1 in subAtoms and a2 in subAtoms:
            numRotBonds+=1
    return numRotBonds


substituentDescriptor = namedtuple('substituentDescriptor',
                                   ['size','relSize','numNO','relNumNO','relNumNO_2','pathLength','relPathLength','relPathLength_2',
                                    'sharedNeighbors', 'numRotBonds', 'numAroBonds'])    

def calcSizeSubstituents(mol, sub_dict, sharedNeighbors_dict, maxShell_dict):  
    sizeDict=defaultdict()
    numAtoms = mol.GetNumAtoms()
    #print(sharedNeighbors_dict)
    for sidx, sub in sorted(sub_dict.items(), key= lambda x: len(x[1])):
        size = _getSizeOfSubstituents(sub, sharedNeighbors_dict)
        numNOs=0
        numShared=0
        # determine the number of oxygen and nitrogen atoms
        for i in sub:
            if mol.GetAtomWithIdx(i).GetAtomicNum() in [7,8]:
                numNOs+=1.0/sharedNeighbors_dict[i]
            if sharedNeighbors_dict[i] > 1:
                numShared+=1
        numRotBs = getRotatableBondsSubstituent(mol, set(sub))
        aroBonds = getAromaticBondsSubstituent(mol, set(sub))
        # fill the substituentDescriptor tuple
        sizeDict[sidx]=substituentDescriptor(size=size,relSize=size/numAtoms,numNO=numNOs,relNumNO=numNOs/numAtoms,
                                             relNumNO_2=numNOs/size,pathLength=maxShell_dict[sidx],
                                             relPathLength=maxShell_dict[sidx]/numAtoms,relPathLength_2=maxShell_dict[sidx]/size,
                                                 sharedNeighbors=numShared, numRotBonds=numRotBs, numAroBonds=aroBonds)
    # if we have less then 4 substituents the missing ones need to be an hydrogen atoms
    if len(sizeDict) < 4:
        for i in range(4-len(sizeDict)):
            sizeDict['H'+str(i)]=substituentDescriptor(size=0,relSize=0,numNO=0,relNumNO=0,relNumNO_2=0,
                                                       pathLength=0,relPathLength=0,relPathLength_2=0,sharedNeighbors=0, numRotBonds=0,
                                                       numAroBonds=0)
    return sizeDict


# Visualization of the substituents
def visualizeSubstituentsGrid(mol, aIdx, molSize=(300,150), kekulize=True,):
    dists = Chem.GetDistanceMatrix(mol)
    idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
    subDict, sn_dict, maxShell_dict = determineAtomSubstituents(aIdx, mol, dists, False)
    
    colors = sns.husl_palette(len(subDict), s=.6)    
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
    count=0
    svgs=[]
    labels=[]
    for sidx, sub in sorted(subDict.items(), key= lambda x: _getSizeOfSubstituents(x[1], sn_dict)):
        color = tuple(colors[count])        
        count+=1
        atColors={}
        for atom in sub:
            atColors[atom]=color
        
        bonds = getBondsSubstituent(mol, set(sub)) 
        bnColors={}
        for b in bonds:
            bnColors[b]=color
        
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
        drawer.DrawMolecule(mc,highlightAtoms=atColors.keys(),
                            highlightAtomColors=atColors,highlightBonds=bonds,highlightBondColors=bnColors)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svgs.append(svg.replace('svg:',''))
        labels.append("Substituent "+str(count)+" (#atoms: "+str(len(sub))+", size normed: "+str(_getSizeOfSubstituents(sub, sn_dict))+")")
    return _svgsToGrid(svgs, labels, svgsPerRow=len(svgs),molSize=molSize,fontSize=12)    
    
def visualizeChiralSubstituentsGrid(mol):
    idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
    return visualizeSubstituentsGrid(mol, idxChiral)

# Chiral moment descriptor    
from numpy import linalg as LA
def calcSP3CarbonSubstituentMoment(subSizes):

    if len(subSizes) != 4:
        raise ValueError('Function "calcSP3CarbonSubstituentMoment" expects an array of size 4 as parameter')
        
    # tetrahedron unit vectors
    x1=np.array([1,1,1])
    x2=np.array([-1,1,-1])
    x3=np.array([1,-1,-1])
    x4=np.array([-1,-1,1])

    substituentMoment= LA.norm((subSizes[0]*x1)+(subSizes[1]*x2)+(subSizes[2]*x3)+(subSizes[3]*x4))
    return substituentMoment

def generateChiralDescriptors(mol):
    desc={}
    dists = Chem.GetDistanceMatrix(mol)
    idxChiral = Chem.FindMolChiralCenters(mol)[0][0]
    subDict, sn_dict, maxShell_dict = determineAtomSubstituents(idxChiral, mol, dists, False)
    sizes = calcSizeSubstituents(mol, subDict, sn_dict, maxShell_dict)
    paths = dists[idxChiral]
    # set some basic descriptors
    desc['numAtoms'] = mol.GetNumAtoms()
    desc['numBonds'] = mol.GetNumBonds()
    desc['numRotBonds'] = AllChem.CalcNumRotatableBonds(mol)
    desc['ringChiralCenter'] = int(mol.GetAtomWithIdx(idxChiral).IsInRing())
    # determine the max path length in the molecule and the mean pairwise distance of all atom pairs
    desc['meanDist'] = np.sum(dists)/((desc['numAtoms']-1)*(desc['numAtoms']))
    desc['maxDist'] = int(np.max(dists))
    # determine the max path length from the chiral center and the mean pairwise distance of 
    # all atom pairs from the chiral center
    desc['meanDistFromCC'] = np.sum(paths)/(desc['numAtoms']-1)
    desc['maxDistfromCC'] = int(np.max(paths))
    # determine the number of neighbors per shell/distance level
    nlevels=Counter(paths.astype(int))
    # consider the levels until a path lenght of 10
    for i in range(1,11):
        desc['nLevel'+str(i)]=nlevels[i]
    # determine the number of nitrogen and oxygen atoms in a certain level around the chiral center
    for i in range(1,4):
        desc['phLevel'+str(i)]=len([n for n,j in enumerate(paths) if j==i and mol.GetAtomWithIdx(n).GetAtomicNum() in [7,8]])
    # determine the number of arometic atoms in a certain level around the chiral center
    for i in range(1,4):
        desc['arLevel'+str(i)]=len([n for n,j in enumerate(paths) if j==i and mol.GetAtomWithIdx(n).GetIsAromatic()])
    # set the size descriptors for each substituent, sort them from smallest to largest
    for n,v in enumerate(sorted(sizes.items(), key=lambda x: x[1].size)):
        desc['s'+str(n+1)+'_size'] = v[1].size
        desc['s'+str(n+1)+'_relSize'] = v[1].relSize
        desc['s'+str(n+1)+'_phSize'] = v[1].numNO
        desc['s'+str(n+1)+'_phRelSize'] = v[1].relNumNO
        desc['s'+str(n+1)+'_phRelSize_2'] = v[1].relNumNO_2
        desc['s'+str(n+1)+'_pathLength'] = v[1].pathLength
        desc['s'+str(n+1)+'_relPathLength'] = v[1].relPathLength
        desc['s'+str(n+1)+'_relPathLength_2'] = v[1].relPathLength_2
        desc['s'+str(n+1)+'_numSharedNeighbors']=v[1].sharedNeighbors
        desc['s'+str(n+1)+'_numRotBonds']=v[1].numRotBonds
        desc['s'+str(n+1)+'_numAroBonds']=v[1].numAroBonds
    # some combination of substituent sizes
    desc['s34_size'] = desc['s3_size']+desc['s4_size']
    desc['s34_phSize'] = desc['s3_phSize']+desc['s4_phSize']
    desc['s34_relSize'] = desc['s3_relSize']+desc['s4_relSize']
    desc['s34_phRelSize'] = desc['s3_phRelSize']+desc['s4_phRelSize']
    # calculate the chiral moment --> kind of 3D descriptor
    desc['chiralMoment'] = calcSP3CarbonSubstituentMoment([desc['s1_size'],desc['s2_size'],desc['s3_size'],desc['s4_size']])
    desc['chiralPhMoment'] = calcSP3CarbonSubstituentMoment([desc['s1_phSize'],desc['s2_phSize'],
                                                                         desc['s3_phSize'],desc['s4_phSize']])
    return desc
