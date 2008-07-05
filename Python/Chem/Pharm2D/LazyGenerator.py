## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" lazy generator of 2D pharmacophore signature data

"""
import Chem
from Chem.Pharm2D import SigFactory,Matcher

class Generator(object):
  """

  Important attributes:

   - mol: the molecules whose signature is being worked with

   - sig: the signature to be used to decide what bits are set in the
     molecule.  This signature can be considered to be constant in
     this class, none of the class methods modify sig

  """
  def __init__(self,sig,mol,dMat=None,bitCache=1):
    """ constructor

      **Arguments**

       - sig: a signature, see class docs

       - mol: a molecule, see class docs

       - dMat: (optional) a distance matrix for the molecule.  If this
         is not provided, one will be calculated

       - bitCache: (optional) if nonzero, a local cache of which bits
         have been queried will be maintained.  Otherwise things must
         be recalculate each time a bit is queried.
       
    """
    if isinstance(sig,SigFactory.SigFactory):
      sig = sig.GetSignature()
    self.sig = sig
    self.mol = mol

    if dMat is None:
      useBO = self.sig.GetIncludeBondOrder()
      dMat = Chem.GetDistanceMatrix(mol,useBO)

    self.dMat = dMat

    if bitCache:
      self.bits = {}
    else:
      self.bits = None

    nPatts = self.sig.GetNumPatterns()
    pattMatches = [None]*nPatts
    for i in range(nPatts):
      patt = self.sig.GetPattern(i)
      pattMatches[i]=self.mol.GetSubstructMatches(patt)
    self.pattMatches = pattMatches
    
  def GetBit(self,idx):
    """ returns a bool indicating whether or not the bit is set

    """
    if idx < 0 or idx >= self.sig.GetSize():
      raise IndexError,'Index %d invalid'%(idx)
    if self.bits is not None and self.bits.has_key(idx):
      return self.bits[idx]
    
    tmp = Matcher.GetAtomsMatchingBit(self.sig,idx,self.mol,
                                      dMat=self.dMat,justOne=1,
                                      matchingAtoms=self.pattMatches)
    if not tmp or len(tmp)==0: res = 0
    else: res = 1
    
    if self.bits is not None:
      self.bits[idx] = res
    return res

  def __len__(self):
    """ allows class to support len()

    """
    return self.sig.GetSize()
  def __getitem__(self,itm):
    """ allows class to support random access.
      Calls self.GetBit()

    """
    return self.GetBit(itm)

    
  

if __name__ == '__main__':
  import time
  import RDConfig,Chem
  from Chem.Pharm2D import Gobbi_Pharm2D,Generate
  import random

  factory = Gobbi_Pharm2D.factory
  nToDo=100
  inD = open(RDConfig.RDDataDir+"/NCI/first_5K.smi",'r').readlines()[:nToDo]
  mols = [None]*len(inD)
  for i in range(len(inD)):
    smi = inD[i].split('\t')[0]
    smi.strip()
    mols[i] = Chem.MolFromSmiles(smi)

  sig = factory.GetSignature()

  nBits = 300
  random.seed(23)
  bits = [random.randint(0,sig.GetSize()-1) for x in range(nBits)]

  print 'Using the Lazy Generator'
  t1 = time.time()
  for i in range(len(mols)):
    if not i % 10: print 'done mol %d of %d'%(i,len(mols))
    gen = Generator(factory,mols[i])
    for bit in bits:
      v = gen[bit]
  t2 = time.time()    
  print '\tthat took %4.2f seconds'%(t2-t1)
    
  
  print 'Generating and checking signatures'
  t1 = time.time()
  for i in range(len(mols)):
    if not i % 10: print 'done mol %d of %d'%(i,len(mols))
    sig = Generate.Gen2DFingerprint(mols[i],factory)
    for bit in bits:
      v = sig[bit]
  t2 = time.time()    
  print '\tthat took %4.2f seconds'%(t2-t1)
  
