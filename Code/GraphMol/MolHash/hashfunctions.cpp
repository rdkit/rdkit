/*==============================================*/
/* Copyright (C)  2011-2019  NextMove Software  */
/* All rights reserved.                         */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#define _CRT_SECURE_NO_WARNINGS

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "nmmolhash.h"
#include "mf.h"

namespace {
unsigned int NMRDKitBondGetOrder(const RDKit::Bond *bnd) {
  PRECONDITION(bnd, "bad bond");
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

RDKit::Bond *NMRDKitMolNewBond(RDKit::RWMol *mol, RDKit::Atom *src,
                               RDKit::Atom *dst, unsigned int order,
                               bool arom) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(src, "bad src atom");
  PRECONDITION(dst, "bad dest atom");
  RDKit::Bond *result;
  result = mol->getBondBetweenAtoms(src->getIdx(), dst->getIdx());
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
      case 1:
        type = RDKit::Bond::SINGLE;
        break;
      case 2:
        type = RDKit::Bond::DOUBLE;
        break;
      case 3:
        type = RDKit::Bond::TRIPLE;
        break;
      case 4:
        type = RDKit::Bond::QUADRUPLE;
        break;
    }
  } else {
    type = RDKit::Bond::AROMATIC;
  }

  result = new RDKit::Bond(type);
  result->setOwningMol(mol);
  result->setBeginAtom(src);
  result->setEndAtom(dst);
  mol->addBond(result, true);
  if (arom) {
    result->setIsAromatic(true);
  }
  return result;
}

void NMRDKitSanitizeHydrogens(RDKit::RWMol *mol) {
  PRECONDITION(mol, "bad molecule");
  // Move all of the implicit Hs into one box
  for (auto aptr : mol->atoms()) {
    unsigned int hcount = aptr->getTotalNumHs();
    aptr->setNoImplicit(true);
    aptr->setNumExplicitHs(hcount);
    aptr->updatePropertyCache();  // or else the valence is reported incorrectly
  }
}

}  // namespace

namespace RDKit {
namespace MolHash {
static unsigned int NMDetermineComponents(RWMol *mol, unsigned int *parts,
                                          unsigned int acount) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(parts, "bad parts pointer");
  memset(parts, 0, acount * sizeof(unsigned int));
  std::vector<Atom *> todo;

  unsigned int result = 0;
  for (auto aptr : mol->atoms()) {
    unsigned int idx = aptr->getIdx();
    if (parts[idx] == 0) {
      parts[idx] = ++result;
      todo.push_back(aptr);

      while (!todo.empty()) {
        aptr = todo.back();
        todo.pop_back();
        for (auto nptr : aptr->nbrs()) {
          idx = nptr->getIdx();
          if (parts[idx] == 0) {
            parts[idx] = result;
            todo.push_back(nptr);
          }
        }
      }
    }
  }
  return result;
}

static std::string NMMolecularFormula(RWMol *mol, const unsigned int *parts,
                                      unsigned int part) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION((!part || parts), "bad parts pointer");
  unsigned int hist[256];
  int charge = 0;

  memset(hist, 0, sizeof(hist));

  for (auto aptr : mol->atoms()) {
    unsigned int idx = aptr->getIdx();
    if (part == 0 || parts[idx] == part) {
      unsigned int elem = aptr->getAtomicNum();
      if (elem < 256) {
        hist[elem]++;
      } else {
        hist[0]++;
      }
      hist[1] += aptr->getTotalNumHs(false);
      charge += aptr->getFormalCharge();
    }
  }

  char buffer[16];
  std::string result;
  const unsigned char *perm = hist[6] ? OrganicHillOrder : InorganicHillOrder;
  for (unsigned int i = 0; i < 119; i++) {
    unsigned int elem = perm[i];
    if (hist[elem] > 0) {
      result += symbol[elem];
      if (hist[elem] > 1) {
        sprintf(buffer, "%u", hist[elem]);
        result += buffer;
      }
    }
  }
  if (charge != 0) {
    if (charge > 0) {
      result += '+';
      if (charge > 1) {
        sprintf(buffer, "%d", charge);
        result += buffer;
      }
    } else {  // charge < 0
      result += '-';
      if (charge < -1) {
        sprintf(buffer, "%d", -charge);
        result += buffer;
      }
    }
  }
  return result;
}

static std::string NMMolecularFormula(RWMol *mol, bool sep = false) {
  PRECONDITION(mol, "bad molecule");
  if (!sep) {
    return NMMolecularFormula(mol, nullptr, 0);
  }

  unsigned int acount = mol->getNumAtoms();
  if (acount == 0) {
    return "";
  }

  auto size = (unsigned int)(acount * sizeof(int));
  auto *parts = (unsigned int *)malloc(size);
  unsigned int pcount = NMDetermineComponents(mol, parts, acount);

  std::string result;
  if (pcount > 1) {
    std::vector<std::string> vmf;
    for (unsigned int i = 1; i <= pcount; i++) {
      vmf.push_back(NMMolecularFormula(mol, parts, i));
    }

    // sort
    result = vmf[0];
    for (unsigned int i = 1; i < pcount; i++) {
      result += ".";
      result += vmf[i];
    }
  } else {  // pcount == 1
    result = NMMolecularFormula(mol, parts, 1);
  }
  free(parts);
  return result;
}

static void NormalizeHCount(Atom *aptr) {
  PRECONDITION(aptr, "bad atom pointer");
  unsigned int hcount;

  switch (aptr->getAtomicNum()) {
    case 9:   // Fluorine
    case 17:  // Chlorine
    case 35:  // Bromine
    case 53:  // Iodine
      hcount = aptr->getDegree();
      hcount = hcount < 1 ? 1 - hcount : 0;
      break;
    case 8:   // Oxygen
    case 16:  // Sulfur
      hcount = aptr->getDegree();
      hcount = hcount < 2 ? 2 - hcount : 0;
      break;
    case 5:  // Boron
      hcount = aptr->getDegree();
      hcount = hcount < 3 ? 3 - hcount : 0;
      break;
    case 7:   // Nitogen
    case 15:  // Phosphorus
      hcount = aptr->getDegree();
      if (hcount < 3) {
        hcount = 3 - hcount;
      } else if (hcount == 4) {
        hcount = 1;
      } else {
        hcount = 0;
      }
      break;
    case 6:  // Carbon
      hcount = aptr->getDegree();
      hcount = hcount < 4 ? 4 - hcount : 0;
      break;
    default:
      hcount = 0;
  }
  aptr->setNoImplicit(true);
  aptr->setNumExplicitHs(hcount);
}

static std::string AnonymousGraph(RWMol *mol, bool elem) {
  PRECONDITION(mol, "bad molecule");
  std::string result;
  int charge = 0;

  for (auto aptr : mol->atoms()) {
    charge += aptr->getFormalCharge();
    aptr->setIsAromatic(false);
    aptr->setFormalCharge(0);
    if (!elem) {
      aptr->setNumExplicitHs(0);
      aptr->setNoImplicit(true);
      aptr->setAtomicNum(0);
    } else {
      NormalizeHCount(aptr);
    }
  }

  for (auto bptr : mol->bonds()) {
    bptr->setBondType(Bond::SINGLE);
  }
  MolOps::assignRadicals(*mol);
  result = MolToSmiles(*mol);
  return result;
}

static std::string MesomerHash(RWMol *mol, bool netq) {
  PRECONDITION(mol, "bad molecule");
  std::string result;
  char buffer[32];
  int charge = 0;

  for (auto aptr : mol->atoms()) {
    charge += aptr->getFormalCharge();
    aptr->setIsAromatic(false);
    aptr->setFormalCharge(0);
  }

  for (auto bptr : mol->bonds()) {
    bptr->setBondType(Bond::SINGLE);
  }

  MolOps::assignRadicals(*mol);
  result = MolToSmiles(*mol);
  if (netq) {
    sprintf(buffer, "_%d", charge);
    result += buffer;
  }
  return result;
}

static std::string TautomerHash(RWMol *mol, bool proto) {
  PRECONDITION(mol, "bad molecule");
  std::string result;
  char buffer[32];
  int hcount = 0;
  int charge = 0;

  for (auto aptr : mol->atoms()) {
    charge += aptr->getFormalCharge();
    aptr->setIsAromatic(false);
    aptr->setFormalCharge(0);
    if (aptr->getAtomicNum() != 6) {
      hcount += aptr->getTotalNumHs(false);
      aptr->setNoImplicit(true);
      aptr->setNumExplicitHs(0);
    }
  }

  for (auto bptr : mol->bonds()) {
    if (bptr->getBondType() != Bond::SINGLE &&
        (bptr->getIsConjugated() || bptr->getBeginAtom()->getAtomicNum() != 6 ||
         bptr->getEndAtom()->getAtomicNum() != 6)) {
      bptr->setIsAromatic(false);
      bptr->setBondType(Bond::SINGLE);
      bptr->setStereo(Bond::BondStereo::STEREONONE);
    }
  }

  MolOps::assignRadicals(*mol);
  // we may have just destroyed some stereocenters/bonds
  // clean that up:
  bool cleanIt = true;
  bool force = true;
  MolOps::assignStereochemistry(*mol, cleanIt, force);
  result = MolToSmiles(*mol);
  if (!proto) {
    sprintf(buffer, "_%d_%d", hcount, charge);
  } else {
    sprintf(buffer, "_%d", hcount - charge);
  }
  result += buffer;
  return result;
}

static bool TraverseForRing(Atom *atom, unsigned char *visit) {
  PRECONDITION(atom, "bad atom pointer");
  PRECONDITION(visit, "bad pointer");
  visit[atom->getIdx()] = 1;
  for (auto nptr : atom->nbrs()) { 
    if (visit[nptr->getIdx()] == 0) {
      if (RDKit::queryIsAtomInRing(nptr)) {
        return true;
      }

      if (TraverseForRing(nptr, visit)) {
        return true;
      }
    }
  }
  return false;
}

static bool DepthFirstSearchForRing(Atom *root, Atom *nbor,
                                    unsigned int maxatomidx) {
  PRECONDITION(root, "bad atom pointer");
  PRECONDITION(nbor, "bad atom pointer");

  unsigned int natoms = maxatomidx;
  auto *visit = (unsigned char *)alloca(natoms);
  memset(visit, 0, natoms);

  visit[root->getIdx()] = true;
  return TraverseForRing(nbor, visit);
}

bool IsInScaffold(Atom *atom, unsigned int maxatomidx) {
  PRECONDITION(atom, "bad atom pointer");
  if (RDKit::queryIsAtomInRing(atom)) {
    return true;
  }

  unsigned int count = 0;
  for (auto nptr : atom->nbrs()) {
    if (DepthFirstSearchForRing(atom, nptr, maxatomidx)) {
      ++count;
    }
  }
  return count > 1;
}

static bool HasNbrInScaffold(Atom *aptr, unsigned char *is_in_scaffold) {
  PRECONDITION(aptr, "bad atom pointer");
  PRECONDITION(is_in_scaffold, "bad pointer");
  for (auto nptr : aptr->nbrs()) {
    if (is_in_scaffold[nptr->getIdx()]) {
      return true;
    }
  }
  return false;
}

static std::string ExtendedMurckoScaffold(RWMol *mol) {
  PRECONDITION(mol, "bad molecule");
  RDKit::MolOps::fastFindRings(*mol);

  unsigned int maxatomidx = mol->getNumAtoms();
  auto *is_in_scaffold = (unsigned char *)alloca(maxatomidx);
  for (auto aptr : mol->atoms()) {
    is_in_scaffold[aptr->getIdx()] = IsInScaffold(aptr, maxatomidx);
  }

  std::vector<Atom *> for_deletion;
  for (auto aptr : mol->atoms()) {
    unsigned int aidx = aptr->getIdx();
    if (is_in_scaffold[aidx]) {
      continue;
    }
    if (HasNbrInScaffold(aptr, is_in_scaffold)) {
      aptr->setAtomicNum(0);
      aptr->setFormalCharge(0);
      aptr->setNoImplicit(true);
      aptr->setNumExplicitHs(0);
    } else {
      for_deletion.push_back(aptr);
    }
  }

  for (auto &i : for_deletion) {
    mol->removeAtom(i);
  }

  MolOps::assignRadicals(*mol);
  std::string result;
  result = MolToSmiles(*mol);
  return result;
}

static std::string MurckoScaffoldHash(RWMol *mol) {
  PRECONDITION(mol, "bad molecule");
  std::vector<Atom *> for_deletion;
  do {
    for_deletion.clear();
    for (auto aptr : mol->atoms()) {
      unsigned int deg = aptr->getDegree();
      if (deg < 2) {
        if (deg == 1) {  // i.e. not 0 and the last atom in the molecule
          for (auto *bptr : aptr->bonds()) {
            Atom *nbr = bptr->getOtherAtom(aptr);
            unsigned int hcount = nbr->getTotalNumHs(false);
            nbr->setNumExplicitHs(hcount + NMRDKitBondGetOrder(bptr));
            nbr->setNoImplicit(true);
          }
        }
        for_deletion.push_back(aptr);
      }
    }
    for (auto &i : for_deletion) {
      mol->removeAtom(i);
    }
  } while (!for_deletion.empty());
  MolOps::assignRadicals(*mol);
  std::string result;
  result = MolToSmiles(*mol);
  return result;
}

static std::string NetChargeHash(RWMol *mol) {
  PRECONDITION(mol, "bad molecule");
  int totalq = 0;

  for (auto aptr : mol->atoms()) {
    totalq += aptr->getFormalCharge();
  }

  char buffer[16];
  sprintf(buffer, "%d", totalq);
  return buffer;
}

static std::string SmallWorldHash(RWMol *mol, bool brl) {
  PRECONDITION(mol, "bad molecule");
  char buffer[64];

  unsigned int acount = mol->getNumAtoms();
  unsigned int bcount = mol->getNumBonds();
  unsigned int rcount = (bcount + 1) - acount;

  if (brl) {
    unsigned int lcount = 0;
    for (auto aptr : mol->atoms()) {
      if (aptr->getDegree() == 2) {
        lcount++;
      }
    }
    sprintf(buffer, "B%uR%uL%u", bcount, rcount, lcount);
  } else {
    sprintf(buffer, "B%uR%u", bcount, rcount);
  }
  return buffer;
}

static void DegreeVector(RWMol *mol, unsigned int *v) {
  memset(v, 0, 4 * sizeof(unsigned int));
  for (auto aptr : mol->atoms()) {
    switch (aptr->getDegree()) {
      case 4:
        v[0]++;
        break;
      case 3:
        v[1]++;
        break;
      case 2:
        v[2]++;
        break;
      case 1:
        v[3]++;
        break;
    }
  }
}

static bool HasDoubleBond(Atom *atom) {
  PRECONDITION(atom, "bad atom");
  for (auto *bptr : atom->bonds()) {
    if (NMRDKitBondGetOrder(bptr) == 2) {
      return true;
    }
  }
  return false;
}

// Determine whether/how to fragment bond
// -1 means don't fragment bond
// 0 means break, with hydrogens on both beg and end
// 1 means break, with asterisk on beg and hydrogen on end
// 2 means break, with hydrogen on beg and asterisk on end
// 3 means break, with asterisks on both beg and end

static int RegioisomerBond(Bond *bnd) {
  PRECONDITION(bnd, "bad bond");
  if (NMRDKitBondGetOrder(bnd) != 1) {
    return -1;
  }
  if (RDKit::queryIsBondInRing(bnd)) {
    return -1;
  }

  Atom *beg = bnd->getBeginAtom();
  Atom *end = bnd->getEndAtom();
  unsigned int beg_elem = beg->getAtomicNum();
  unsigned int end_elem = end->getAtomicNum();

  if (beg_elem == 0 || end_elem == 0) {
    return -1;
  }

  if (RDKit::queryIsAtomInRing(beg)) {
    if (RDKit::queryIsAtomInRing(end)) {
      return 0;
    }
    return 2;
  }
  if (RDKit::queryIsAtomInRing(end)) {
    return 1;
  }

  if (beg_elem != 6 && end_elem == 6 && !HasDoubleBond(end)) {
    return 1;
  }
  if (beg_elem == 6 && end_elem != 6 && !HasDoubleBond(beg)) {
    return 2;
  }

  return -1;
}

static void ClearEZStereo(Atom *atm) {
  PRECONDITION(atm, "bad atom");
  for (auto *bptr : atm->bonds()) {
    if (bptr->getStereo() > RDKit::Bond::STEREOANY) {
      bptr->setStereo(RDKit::Bond::STEREOANY);
    }
  }
}

static std::string RegioisomerHash(RWMol *mol) {
  PRECONDITION(mol, "bad molecule");

  // we need a copy of the molecule so that we can loop over the bonds of
  // something while modifying something else
  RDKit::MolOps::fastFindRings(*mol);
  RDKit::ROMol molcpy(*mol);
  for (int i = molcpy.getNumBonds() - 1; i >= 0; --i) {
    auto bptr = molcpy.getBondWithIdx(i);
    int split = RegioisomerBond(bptr);
    if (split >= 0) {
      bptr = mol->getBondWithIdx(i);
      Atom *beg = bptr->getBeginAtom();
      Atom *end = bptr->getEndAtom();
      mol->removeBond(bptr->getBeginAtomIdx(), bptr->getEndAtomIdx());
      ClearEZStereo(beg);
      ClearEZStereo(end);

      if (split & 1) {
        Atom *star = new RDKit::Atom(0);
        mol->addAtom(star, true, true);
        star->setNoImplicit(true);
        NMRDKitMolNewBond(mol, beg, star, 1, false);
      } else {
        unsigned int hcount = beg->getTotalNumHs(false);
        beg->setNumExplicitHs(hcount + 1);
        beg->setNoImplicit(true);
      }
      if (split & 2) {
        Atom *star = new RDKit::Atom(0);
        mol->addAtom(star, true, true);
        star->setNoImplicit(true);
        NMRDKitMolNewBond(mol, end, star, 1, false);
      } else {
        unsigned int hcount = end->getTotalNumHs(false);
        end->setNumExplicitHs(hcount + 1);
        end->setNoImplicit(true);
      }
    }
  }

  std::string result;
  result = MolToSmiles(*mol);
  return result;
}

static std::string ArthorSubOrderHash(RWMol *mol) {
  PRECONDITION(mol, "bad molecule");
  char buffer[256];

  unsigned int acount = mol->getNumAtoms();
  unsigned int bcount = mol->getNumBonds();

  unsigned int pcount = 1;
  unsigned int size = 4 * mol->getNumAtoms() + 4;
  auto *parts = (unsigned int *)malloc(size);
  if (parts) {
    memset(parts, 0, size);
    pcount = NMDetermineComponents(mol, parts, acount);
    free(parts);
  }

  unsigned int ccount = 0;
  unsigned int ocount = 0;
  unsigned int zcount = 0;
  unsigned int icount = 0;
  unsigned int qcount = 0;
  unsigned int rcount = 0;

  for (auto aptr : mol->atoms()) {
    unsigned int elem = aptr->getAtomicNum();
    int charge = aptr->getFormalCharge();
    switch (elem) {
      case 6:  // Carbon
        ccount++;
        if (charge == 0 && aptr->getTotalValence() != 4) {
          rcount++;
        }
        break;
      case 7:   // Nitrogen
      case 15:  // Phosphorus
        ocount++;
        if (charge == 0) {
          unsigned int valence = aptr->getTotalValence();
          if (valence != 3 && valence != 5) {
            rcount++;
          }
        }
        break;
      case 8:  // Oxygen
        ocount++;
        if (charge && aptr->getTotalValence() != 2) {
          rcount++;
        }
        break;
      case 9:  // Fluorine
        ocount++;
        if (charge && aptr->getTotalValence() != 1) {
          rcount++;
        }
        break;
      case 17:  // Chlorine
      case 35:  // Bromine
      case 53:  // Iodine
        ocount++;
        if (charge == 0) {
          unsigned int valence = aptr->getTotalValence();
          if (valence != 1 && valence != 3 && valence != 5 && valence != 7) {
            rcount++;
          }
        }
        break;
      case 16:  // Sulfur
        ocount++;
        if (charge == 0) {
          unsigned int valence = aptr->getTotalValence();
          if (valence != 2 && valence != 4 && valence != 6) {
            rcount++;
          }
        }
        break;
    }
    zcount += elem;
    if (aptr->getIsotope() != 0) {
      icount++;
    }
    if (charge != 0) {
      qcount++;
    }
  }

  if (acount > 0xffff) {
    acount = 0xffff;
  }
  if (bcount > 0xffff) {
    bcount = 0xffff;
  }
  if (pcount > 0xff) {
    pcount = 0xff;
  }
  if (ccount > 0xffff) {
    ccount = 0xffff;
  }
  if (ocount > 0xffff) {
    ocount = 0xffff;
  }
  if (zcount > 0xffffff) {
    zcount = 0xffffff;
  }
  if (rcount > 0xff) {
    rcount = 0xff;
  }
  if (qcount > 0xff) {
    qcount = 0xff;
  }
  if (icount > 0xff) {
    icount = 0xff;
  }

  sprintf(buffer, "%04x%04x%02x%04x%04x%06x%02x%02x%02x", acount, bcount,
          pcount, ccount, ocount, zcount, rcount, qcount, icount);
  return buffer;
}

std::string MolHash(RWMol *mol, HashFunction func) {
  PRECONDITION(mol, "bad molecule");
  std::string result;
  char buffer[32];
  NMRDKitSanitizeHydrogens(mol);

  switch (func) {
    default:
    case HashFunction::AnonymousGraph:
      result = AnonymousGraph(mol, false);
      break;
    case HashFunction::ElementGraph:
      result = AnonymousGraph(mol, true);
      break;
    case HashFunction::CanonicalSmiles:
      result = MolToSmiles(*mol);
      break;
    case HashFunction::MurckoScaffold:
      result = MurckoScaffoldHash(mol);
      break;
    case HashFunction::ExtendedMurcko:
      result = ExtendedMurckoScaffold(mol);
      break;
    case HashFunction::Mesomer:
      result = MesomerHash(mol, true);
      break;
    case HashFunction::RedoxPair:
      result = MesomerHash(mol, false);
      break;
    case HashFunction::HetAtomTautomer:
      result = TautomerHash(mol, false);
      break;
    case HashFunction::HetAtomProtomer:
      result = TautomerHash(mol, true);
      break;
    case HashFunction::MolFormula:
      result = NMMolecularFormula(mol);
      break;
    case HashFunction::AtomBondCounts:
      sprintf(buffer, "%u,%u", mol->getNumAtoms(), mol->getNumBonds());
      result = buffer;
      break;
    case HashFunction::NetCharge:
      result = NetChargeHash(mol);
      break;
    case HashFunction::SmallWorldIndexBR:
      result = SmallWorldHash(mol, false);
      break;
    case HashFunction::SmallWorldIndexBRL:
      result = SmallWorldHash(mol, true);
      break;
    case HashFunction::DegreeVector: {
      unsigned int dv[4];
      DegreeVector(mol, dv);
      sprintf(buffer, "%u,%u,%u,%u", dv[0], dv[1], dv[2], dv[3]);
      result = buffer;
    } break;
    case HashFunction::ArthorSubstructureOrder:
      result = ArthorSubOrderHash(mol);
      break;
    case HashFunction::Regioisomer:
      result = RegioisomerHash(mol);
      break;
  }
  return result;
}
}  // namespace MolHash
}  // namespace RDKit
