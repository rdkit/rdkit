//
//  Copyright (C) 2013-2017 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ProximityBonds.h"
#include <algorithm>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MonomerInfo.h>

namespace RDKit {

static const double EXTDIST = 0.45;
// static const double MAXRAD   = 2.50;
// static const double MINDIST  = 0.40;
static const double MAXDIST = 5.45;      // 2*MAXRAD + EXTDIST
static const double MINDIST2 = 0.16;     // MINDIST*MINDIST
static const double MAXDIST2 = 29.7025;  // MAXDIST*MAXDIST

struct ProximityEntry {
  float x, y, z, r;
  int atm, hash, next, elem;

  bool operator<(const ProximityEntry &p) const { return x < p.x; }
};

static bool IsBonded(ProximityEntry *p, ProximityEntry *q, unsigned int flags) {
  if (flags & ctdIGNORE_H_H_CONTACTS && p->elem == 1 && q->elem == 1) {
    return false;
  }
  double dx = (double)p->x - (double)q->x;
  double dist2 = dx * dx;
  if (dist2 > MAXDIST2) {
    return false;
  }
  double dy = (double)p->y - (double)q->y;
  dist2 += dy * dy;
  if (dist2 > MAXDIST2) {
    return false;
  }
  double dz = (double)p->z - (double)q->z;
  dist2 += dz * dz;

  if (dist2 > MAXDIST2 || dist2 < MINDIST2) {
    return false;
  }

  double radius = (double)p->r + (double)q->r + EXTDIST;
  return dist2 <= radius * radius;
}

bool SamePDBResidue(AtomPDBResidueInfo *p, AtomPDBResidueInfo *q) {
  return p->getResidueNumber() == q->getResidueNumber() &&
         p->getResidueName() == q->getResidueName() &&
         p->getChainId() == q->getChainId() &&
         p->getInsertionCode() == q->getInsertionCode();
}

static bool IsBlacklistedAtom(Atom *atom) {
  // blacklist metals, noble gasses and halogens
  int elem = atom->getAtomicNum();
  // make an inverse query (non-metals and metaloids)
  return !((5 <= elem && elem <= 8) || (14 <= elem && elem <= 16) ||
           (32 <= elem && elem <= 34) || (51 <= elem && elem <= 52));
}

bool IsBlacklistedPair(Atom *beg_atom, Atom *end_atom) {
  PRECONDITION(beg_atom, "empty atom");
  PRECONDITION(end_atom, "empty atom");

  auto *beg_info = (AtomPDBResidueInfo *)beg_atom->getMonomerInfo();
  auto *end_info = (AtomPDBResidueInfo *)end_atom->getMonomerInfo();
  if (!beg_info || beg_info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return false;
  }
  if (!end_info || end_info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return false;
  }

  if (!SamePDBResidue(beg_info, end_info)) {
    if (IsBlacklistedAtom(beg_atom) || IsBlacklistedAtom(end_atom)) {
      return true;
    }
    // Dont make bonds to waters
    if (beg_info->getResidueName() == "HOH" ||
        end_info->getResidueName() == "HOH") {
      return true;
    }
  }
  return false;
}

/*
static void ConnectTheDots_Small(RWMol *mol)
{
unsigned int count = mol->getNumAtoms();
ProximityEntry *tmp = (ProximityEntry*)malloc(count*sizeof(ProximityEntry));
PeriodicTable *table = PeriodicTable::getTable();
Conformer *conf = &mol->getConformer();
for (unsigned int i=0; i<count; i++) {
  Atom *atom = mol->getAtomWithIdx(i);
  unsigned int elem = atom->getAtomicNum();
  RDGeom::Point3D p = conf->getAtomPos(i);
  ProximityEntry *tmpi = tmp+i;
  tmpi->x = (float)p.x;
  tmpi->y = (float)p.y;
  tmpi->z = (float)p.z;
  tmpi->r = (float)table->getRcovalent(elem);
  for (unsigned int j=0; j<i; j++) {
    ProximityEntry *tmpj = tmp+j;
    if (IsBonded(tmpi,tmpj) && !mol->getBondBetweenAtoms(i,j))
      mol->addBond(i,j,Bond::SINGLE);
  }
}
free(tmp);
}


static void ConnectTheDots_Medium(RWMol *mol)
{
int count = mol->getNumAtoms();
std::vector<ProximityEntry> tmp(count);
PeriodicTable *table = PeriodicTable::getTable();
Conformer *conf = &mol->getConformer();
for (int i=0; i<count; i++) {
  Atom *atom = mol->getAtomWithIdx(i);
  unsigned int elem = atom->getAtomicNum();
  RDGeom::Point3D p = conf->getAtomPos(i);
  ProximityEntry *tmpi = &tmp[i];
  tmpi->x = (float)p.x;
  tmpi->y = (float)p.y;
  tmpi->z = (float)p.z;
  tmpi->r = (float)table->getRcovalent(elem);
  tmpi->atm = i;
}

std::stable_sort(tmp.begin(),tmp.end());

for (int j=0; j<count; j++) {
  ProximityEntry *tmpj = &tmp[j];
  double limit = tmpj->x - MAXDIST;
  for (int k=j-1; k>=0; k--) {
    ProximityEntry *tmpk = &tmp[k];
    if (tmpk->x < limit)
      break;
    if (IsBonded(tmpj,tmpk) &&
        !mol->getBondBetweenAtoms(tmpj->atm,tmpk->atm))
      mol->addBond(tmpj->atm,tmpk->atm,Bond::SINGLE);
  }
}
}
*/

constexpr int HASHSIZE = 1024;
constexpr int HASHMASK = 1023;
constexpr int HASHX = 571;
constexpr int HASHY = 127;
constexpr int HASHZ = 3;

static void ConnectTheDots_Large(RWMol *mol, unsigned int flags) {
  int HashTable[HASHSIZE];
  memset(HashTable, -1, sizeof(HashTable));

  unsigned int count = mol->getNumAtoms();
  auto *tmp = (ProximityEntry *)malloc(count * sizeof(ProximityEntry));
  CHECK_INVARIANT(tmp, "bad allocation");
  PeriodicTable *table = PeriodicTable::getTable();
  Conformer *conf = &mol->getConformer();

  for (unsigned int i = 0; i < count; i++) {
    Atom *atom = mol->getAtomWithIdx(i);
    unsigned int elem = atom->getAtomicNum();
    RDGeom::Point3D p = conf->getAtomPos(i);
    ProximityEntry *tmpi = tmp + i;
    tmpi->x = (float)p.x;
    tmpi->y = (float)p.y;
    tmpi->z = (float)p.z;
    tmpi->r = (float)table->getRcovalent(elem);
    tmpi->atm = i;
    tmpi->elem = elem;

    int hash = HASHX * (int)(p.x / MAXDIST) + HASHY * (int)(p.y / MAXDIST) +
               HASHZ * (int)(p.z / MAXDIST);

    for (int dx = -HASHX; dx <= HASHX; dx += HASHX) {
      for (int dy = -HASHY; dy <= HASHY; dy += HASHY) {
        for (int dz = -HASHZ; dz <= HASHZ; dz += HASHZ) {
          int probe = hash + dx + dy + dz;
          int list = HashTable[probe & HASHMASK];
          while (list != -1) {
            ProximityEntry *tmpj = &tmp[list];
            if (tmpj->hash == probe && IsBonded(tmpi, tmpj, flags) &&
                !mol->getBondBetweenAtoms(tmpi->atm, tmpj->atm) &&
                !IsBlacklistedPair(atom, mol->getAtomWithIdx(tmpj->atm))) {
              mol->addBond(tmpi->atm, tmpj->atm, Bond::SINGLE);
            }
            list = tmpj->next;
          }
        }
      }
    }
    int list = hash & HASHMASK;
    tmpi->next = HashTable[list];
    HashTable[list] = i;
    tmpi->hash = hash;
  }
  // Cleanup pass
  for (unsigned int i = 0; i < count; i++) {
    Atom *atom = mol->getAtomWithIdx(i);
    unsigned int elem = atom->getAtomicNum();
    // detect multivalent Hs, which could happen with ConnectTheDots
    if (elem == 1 && atom->getDegree() > 1) {
      auto *atom_info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
      // cut all but shortest Bond
      RDGeom::Point3D p = conf->getAtomPos(i);
      RDKit::RWMol::ADJ_ITER nbr, end_nbr;
      boost::tie(nbr, end_nbr) = mol->getAtomNeighbors(atom);
      float best = 10000;
      unsigned int best_idx = mol->getNumAtoms() + 1;
      while (nbr != end_nbr) {
        RDGeom::Point3D pn = conf->getAtomPos(*nbr);
        float d = (p - pn).length();
        auto *n_info =
            (AtomPDBResidueInfo *)(mol->getAtomWithIdx(*nbr)->getMonomerInfo());
        if (d < best &&
            atom_info->getResidueNumber() == n_info->getResidueNumber()) {
          best = d;
          best_idx = *nbr;
        }
        ++nbr;
      }
      // iterate again and remove all but closest
      boost::tie(nbr, end_nbr) = mol->getAtomNeighbors(atom);
      while (nbr != end_nbr) {
        if (*nbr == best_idx) {
          Bond *bond = mol->getBondBetweenAtoms(i, *nbr);
          bond->setBondType(Bond::SINGLE);  // make sure this one is single
        } else {
          mol->removeBond(i, *nbr);
        }
        ++nbr;
      }
    }
  }
  free(tmp);
}

void ConnectTheDots(RWMol *mol, unsigned int flags) {
  if (!mol || !mol->getNumConformers()) {
    return;
  }
  // Determine optimal algorithm to use by getNumAtoms()?
  ConnectTheDots_Large(mol, flags);
}

// These are macros to allow their use in C++ constants
constexpr unsigned BCNAM(char A, char B, char C) {
  return (A << 16) | (B << 8) | C;
}
constexpr unsigned BCATM(char A, char B, char C, char D) {
  return (A << 24) | (B << 16) | (C << 8) | D;
}

static bool StandardPDBDoubleBond(unsigned int rescode, unsigned int atm1,
                                  unsigned int atm2) {
  if (atm1 > atm2) {
    unsigned int tmp = atm1;
    atm1 = atm2;
    atm2 = tmp;
  }

  switch (rescode) {
    case BCNAM('A', 'C', 'E'):
    case BCNAM('A', 'L', 'A'):
    case BCNAM('C', 'Y', 'S'):
    case BCNAM('G', 'L', 'Y'):
    case BCNAM('I', 'L', 'E'):
    case BCNAM('L', 'E', 'U'):
    case BCNAM('L', 'Y', 'S'):
    case BCNAM('M', 'E', 'T'):
    case BCNAM('P', 'R', 'O'):
    case BCNAM('S', 'E', 'R'):
    case BCNAM('T', 'H', 'R'):
    case BCNAM('V', 'A', 'L'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      break;
    case BCNAM('A', 'R', 'G'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'Z', ' ') &&
          atm2 == BCATM(' ', 'N', 'H', '2')) {
        return true;
      }
      break;
    case BCNAM('A', 'S', 'N'):
    case BCNAM('A', 'S', 'P'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'G', ' ') &&
          atm2 == BCATM(' ', 'O', 'D', '1')) {
        return true;
      }
      break;
    case BCNAM('G', 'L', 'N'):
    case BCNAM('G', 'L', 'U'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', ' ') &&
          atm2 == BCATM(' ', 'O', 'E', '1')) {
        return true;
      }
      break;
    case BCNAM('H', 'I', 'S'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', '2') &&
          atm2 == BCATM(' ', 'C', 'G', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'E', '1') &&
          atm2 == BCATM(' ', 'N', 'D', '1')) {
        return true;
      }
      break;
    case BCNAM('P', 'H', 'E'):
    case BCNAM('T', 'Y', 'R'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', '1') &&
          atm2 == BCATM(' ', 'C', 'G', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', '2') &&
          atm2 == BCATM(' ', 'C', 'E', '2')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'E', '1') &&
          atm2 == BCATM(' ', 'C', 'Z', ' ')) {
        return true;
      }
      break;
    case BCNAM('T', 'R', 'P'):
      if (atm1 == BCATM(' ', 'C', ' ', ' ') &&
          atm2 == BCATM(' ', 'O', ' ', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', '1') &&
          atm2 == BCATM(' ', 'C', 'G', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'D', '2') &&
          atm2 == BCATM(' ', 'C', 'E', '2')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'E', '3') &&
          atm2 == BCATM(' ', 'C', 'Z', '3')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', 'H', '2') &&
          atm2 == BCATM(' ', 'C', 'Z', '2')) {
        return true;
      }
      break;
    case BCNAM(' ', ' ', 'A'):
    case BCNAM(' ', 'D', 'A'):
      if (atm1 == BCATM(' ', 'C', '6', ' ') &&
          atm2 == BCATM(' ', 'N', '1', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '2', ' ') &&
          atm2 == BCATM(' ', 'N', '3', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '4', ' ') &&
          atm2 == BCATM(' ', 'C', '5', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '8', ' ') &&
          atm2 == BCATM(' ', 'N', '7', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'O', 'P', '1') &&
          atm2 == BCATM(' ', 'P', ' ', ' ')) {
        return true;
      }
      break;
    case BCNAM(' ', ' ', 'G'):
    case BCNAM(' ', 'D', 'G'):
      if (atm1 == BCATM(' ', 'C', '6', ' ') &&
          atm2 == BCATM(' ', 'O', '6', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '2', ' ') &&
          atm2 == BCATM(' ', 'N', '3', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '4', ' ') &&
          atm2 == BCATM(' ', 'C', '5', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '8', ' ') &&
          atm2 == BCATM(' ', 'N', '7', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'O', 'P', '1') &&
          atm2 == BCATM(' ', 'P', ' ', ' ')) {
        return true;
      }
      break;
    case BCNAM(' ', ' ', 'C'):
    case BCNAM(' ', 'D', 'C'):
      if (atm1 == BCATM(' ', 'C', '2', ' ') &&
          atm2 == BCATM(' ', 'O', '2', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '4', ' ') &&
          atm2 == BCATM(' ', 'N', '3', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '5', ' ') &&
          atm2 == BCATM(' ', 'C', '6', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'O', 'P', '1') &&
          atm2 == BCATM(' ', 'P', ' ', ' ')) {
        return true;
      }
      break;
    case BCNAM(' ', ' ', 'T'):
    case BCNAM(' ', 'D', 'T'):
    case BCNAM(' ', ' ', 'U'):
    case BCNAM(' ', 'D', 'U'):
      if (atm1 == BCATM(' ', 'C', '2', ' ') &&
          atm2 == BCATM(' ', 'O', '2', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '4', ' ') &&
          atm2 == BCATM(' ', 'O', '4', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'C', '5', ' ') &&
          atm2 == BCATM(' ', 'C', '6', ' ')) {
        return true;
      }
      if (atm1 == BCATM(' ', 'O', 'P', '1') &&
          atm2 == BCATM(' ', 'P', ' ', ' ')) {
        return true;
      }
      break;
  }
  return false;
}

static bool StandardPDBDoubleBond(RWMol *mol, Atom *beg, Atom *end) {
  auto *bInfo = (AtomPDBResidueInfo *)beg->getMonomerInfo();
  if (!bInfo || bInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return false;
  }
  auto *eInfo = (AtomPDBResidueInfo *)end->getMonomerInfo();
  if (!eInfo || eInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return false;
  }
  if (!SamePDBResidue(bInfo, eInfo)) {
    return false;
  }
  if (bInfo->getIsHeteroAtom() || eInfo->getIsHeteroAtom()) {
    return false;
  }

  const char *ptr = bInfo->getResidueName().c_str();
  unsigned int rescode = BCNAM(ptr[0], ptr[1], ptr[2]);
  ptr = bInfo->getName().c_str();
  unsigned int atm1 = BCATM(ptr[0], ptr[1], ptr[2], ptr[3]);
  ptr = eInfo->getName().c_str();
  unsigned int atm2 = BCATM(ptr[0], ptr[1], ptr[2], ptr[3]);

  if (!StandardPDBDoubleBond(rescode, atm1, atm2)) {
    return false;
  }

  // Check that neither end already has a double bond
  ROMol::OBOND_ITER_PAIR bp;
  for (bp = mol->getAtomBonds(beg); bp.first != bp.second; ++bp.first) {
    if ((*mol)[*bp.first]->getBondType() == Bond::DOUBLE) {
      return false;
    }
  }
  for (bp = mol->getAtomBonds(end); bp.first != bp.second; ++bp.first) {
    if ((*mol)[*bp.first]->getBondType() == Bond::DOUBLE) {
      return false;
    }
  }
  return true;
}

void StandardPDBResidueBondOrders(RWMol *mol) {
  RWMol::BondIterator bondIt;
  for (bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    Bond *bond = *bondIt;
    if (bond->getBondType() == Bond::SINGLE) {
      Atom *beg = bond->getBeginAtom();
      Atom *end = bond->getEndAtom();
      if (StandardPDBDoubleBond(mol, beg, end)) {
        bond->setBondType(Bond::DOUBLE);
      }
    }
  }
}

}  // namespace RDKit
