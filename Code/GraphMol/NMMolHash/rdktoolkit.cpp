/*==============================================*/
/* Copyright (C)  2012-2019  NextMove Software  */
/* All rights reserved.                         */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#include "toolkit.h"

#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>


#include "RDGeneral/versions.h"

RDKit::RWMol *NMRDKitSmilesToMol(const char *smi)
{
  // Special case empty molecules for RDKit
  if (smi[0] == ' ' || smi[0] == '\t' || smi[0] == '\0')
    return new RDKit::RWMol();
  RDKit::RWMol *mol = (RDKit::RWMol*)0;
  try {
    mol = RDKit::SmilesToMol(smi);
  }
  catch (...) { // RDKit::MolSanitizeException and friends
    mol = (RDKit::RWMol*)0;
  }

  // Extract title
  const char* ptr = smi;
  while(true) {
    if (*ptr == '\0' || *ptr == ' ' || *ptr == '\t')
      break;
    ptr++;
  }
  while (*ptr == ' ' || *ptr == '\t')
    ptr++;
  if (*ptr) {
    std::string title(ptr);
    mol->setProp("_Name", title);
  }

  return mol;
}


RDKit::Atom *NMRDKitMolNewAtom(RDKit::RWMol *mol, unsigned int elem)
{
  RDKit::Atom *result = new RDKit::Atom(elem);
  mol->addAtom(result,true,true);
  result->setNoImplicit(true);

#if 0
  /* This should now be fixed in latest RDKit */
  if (mol->getNumConformers())
    mol->getConformer().resize(mol->getNumAtoms());
#endif

  return result;
}


RDKit::Bond *NMRDKitMolNewBond(RDKit::RWMol *mol,
                               RDKit::Atom *src, RDKit::Atom *dst,
                               unsigned int order, bool arom)
{
  RDKit::Bond *result;
  result = mol->getBondBetweenAtoms(src->getIdx(),dst->getIdx());
  if (result) {
    if (order == 1) {
      switch (result->getBondType()) {
      case RDKit::Bond::SINGLE:
        result->setBondType(RDKit::Bond::DOUBLE);
        break;
      case RDKit::Bond::DOUBLE:
        result->setBondType(RDKit::Bond::TRIPLE);
        break;
      default:
        break;
      }
    }
    return result;
  }
  RDKit::Bond::BondType type = RDKit::Bond::UNSPECIFIED;
  if (!arom) {
    switch (order) {
    case 1:  type = RDKit::Bond::SINGLE;     break;
    case 2:  type = RDKit::Bond::DOUBLE;     break;
    case 3:  type = RDKit::Bond::TRIPLE;     break;
    case 4:  type = RDKit::Bond::QUADRUPLE;  break;
    }
  } else type = RDKit::Bond::AROMATIC;

  result = new RDKit::Bond(type);
  result->setOwningMol(mol);
  result->setBeginAtom(src);
  result->setEndAtom(dst);
  mol->addBond(result,true);
  if (arom)
    result->setIsAromatic(true);
  return result;
}


unsigned int NMRDKitAtomGetExplicitValence(RDKit::Atom *atm)
{
  unsigned int result = 0;
  RDKit::ROMol::OEDGE_ITER beg,end;
  const RDKit::ROMol &mol = atm->getOwningMol();
  boost::tie(beg,end) = mol.getAtomBonds(atm);
  while (beg!=end) {
    switch (mol[*beg]->getBondType()) {
    default:  /* silence warnings */
    case RDKit::Bond::SINGLE:
    case RDKit::Bond::AROMATIC:
      result += 1;
      break;
    case RDKit::Bond::DOUBLE:
      result += 2;
      break;
    case RDKit::Bond::TRIPLE:
      result += 3;
      break;
    case RDKit::Bond::QUADRUPLE:
      result += 4;
      break;
    }
    ++beg;
  }
  return result;
}


void NMRDKitBondSetOrder(RDKit::Bond *bnd, unsigned int order)
{
  switch (order) {
  case 1:  bnd->setBondType(RDKit::Bond::SINGLE);  break;
  case 2:  bnd->setBondType(RDKit::Bond::DOUBLE);  break;
  case 3:  bnd->setBondType(RDKit::Bond::TRIPLE);  break;
  case 4:  bnd->setBondType(RDKit::Bond::QUADRUPLE);  break;
  case 5:  bnd->setBondType(RDKit::Bond::QUINTUPLE);  break;
  }
  bnd->setIsAromatic(false);
}


unsigned int NMRDKitBondGetOrder(const RDKit::Bond *bnd)
{
  switch (bnd->getBondType()) {
  case RDKit::Bond::AROMATIC:
  case RDKit::Bond::SINGLE:
    return 1;
  case RDKit::Bond::DOUBLE:
    return 2;
  case RDKit::Bond::TRIPLE:
    return 3;
  case RDKit::Bond::QUADRUPLE:
    return 4;
  case RDKit::Bond::QUINTUPLE:
    return 5;
  case RDKit::Bond::HEXTUPLE:
    return 6;
  default:
    return 0;
  }
}


std::string NMRDKitMolGetTitle(RDKit::RWMol *mol)
{
  static std::string key("_Name");
  std::string result;
  if (mol->hasProp(key))
    mol->getProp(key,result);
  return result;
}


void NMRDKitAtomSetMapIdx(RDKit::Atom *atm, unsigned int idx)
{
  static std::string key("molAtomMapNumber");
  if (idx) 
    atm->setProp(key,(int)idx);
  else if (atm->hasProp(key))
    atm->clearProp(key);
}


void NMRDKitAtomSetImplicitHCount(RDKit::Atom *atm, unsigned int hcount)
{
  atm->setNoImplicit(true);
  atm->setNumExplicitHs(hcount);
}

void NMRDKitSanitizeHydrogens(RDKit::RWMol *mol)
{
  // Move all of the implicit Hs into one box
  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    unsigned int hcount = aptr->getTotalNumHs();
    aptr->setNoImplicit(true);
    aptr->setNumExplicitHs(hcount);
    aptr->updatePropertyCache(); // or else the valence is reported incorrectly
  }
}

void NMRDKitMolSplitFragments(NMS_pMOL mol, std::vector<NMS_MOL*> &fragments)
{
  RDKit::MOL_SPTR_VECT mfrags = RDKit::MolOps::getMolFrags(*mol);
  RDKit::MOL_SPTR_VECT::iterator vit;
  for(vit = mfrags.begin(); vit != mfrags.end(); ++vit) {
    RDKit::ROMol* wrappedmol = (*vit).get(); // reach inside the shared pointer...
    fragments.push_back(new NMS_MOL(*wrappedmol)); // ...and make a copy
  }
}

void NMRDKitMolCalculateRingInfo(NMS_pMOL mol)
{
  if (!mol->getRingInfo()->isInitialized())
    RDKit::MolOps::fastFindRings(*mol);
}
