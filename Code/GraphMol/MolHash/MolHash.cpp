// $Id$
//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <stdio.h>
#include <stdlib.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <boost/crc.hpp>
#include <cstdint>
#include "../Descriptors/MolDescriptors.h"

#include "MolHash.h"

namespace RDKit {
namespace MolHash {
struct MolFragment  // Reference to a fragment of source molecule
    {
  std::vector<const Atom *> Atoms;
  std::vector<const Bond *> Bonds;
  std::vector<std::uint32_t> AtomsIdx;
  std::vector<std::uint32_t> BondsIdx;
  std::map<std::uint32_t, std::uint32_t> MolAtomIdxMap;  // Full Molecule to
                                                             // fragment indeces
                                                             // backward
                                                             // conversion map
 public:
  std::uint32_t getNumAtoms() const { return AtomsIdx.size(); }
  std::uint32_t getNumBonds() const { return BondsIdx.size(); }
};

// INTERNAL FUNCTIONS:
static HashCodeType computeMorganCodeHash(
    const MolFragment &mol, const std::vector<std::uint32_t> &atomLabels,
    const std::vector<std::uint32_t> &bondLabels);
static void prepareMolFragment(MolFragment &m, const ROMol &mol,
                               const std::vector<unsigned> *atomsToUse,
                               const std::vector<unsigned> *bondsToUse);
static void prepareLabels(std::vector<std::uint32_t> &atomLabels,
                          std::vector<std::uint32_t> &bondLabels,
                          const ROMol &mol, const MolFragment &m,
                          const std::vector<std::uint32_t> *atomCodes,
                          const std::vector<std::uint32_t> *bondCodes);
static std::uint32_t computeCRC32(const void *data, size_t size) {
  boost::crc_32_type crc;
  crc.process_bytes(data, size);
  return crc.checksum();
}
//=============================================================================
// MolHash Module API implementation:
//=============================================================================

void fillAtomBondCodes(
    const ROMol &mol, boost::uint64_t flags  // CodeFlags constants combination
    ,
    std::vector<std::uint32_t> *atomCodes  // NULL is allowed
    ,
    std::vector<std::uint32_t> *bondCodes)  // NULL is allowed
{
  if (atomCodes) {
    unsigned n = mol.getNumAtoms();
    atomCodes->resize(n);
    for (unsigned i = 0; i < n; i++) {
      if (0 == (CF_ATOM_ALL & flags))  // NO LABELS
      {
        (*atomCodes)[i] = 1;
        continue;
      }
      const Atom *atom = mol.getAtomWithIdx(i);
      (*atomCodes)[i] = 0;
      if (0 != (CF_ELEMENT & flags)) (*atomCodes)[i] |= atom->getAtomicNum();
      if (0 != (CF_CHARGE & flags))
        (*atomCodes)[i] |= (atom->getFormalCharge() + 8)
                           << 8;  // allowed range [-8, +8]
      if (0 != (CF_VALENCE & flags))
        (*atomCodes)[i] |= (atom->getExplicitValence())
                           << 13;  // getTotalValence()
      if (0 != (CF_ATOM_CHIRALITY & flags)) {
        char v = 0;
        if (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
            atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
          if (atom->hasProp("_CIPCode")) {
            std::string code = atom->getProp<std::string>("_CIPCode");
            if (code == "R")
              v = 1;
            else if (code == "S")
              v = 2;
          } else if (atom->hasProp("_ringStereoAtoms")) {
            const INT_VECT &ringStereoAtoms =
                atom->getProp<INT_VECT>("_ringStereoAtoms");
            if (ringStereoAtoms.size()) {
              if (ringStereoAtoms[0] < 0) {
                v = 1;
              } else {
                v = 2;
              }
              if (ringStereoAtoms.size() > 1) {
                BOOST_LOG(rdWarningLog)
                    << "Warning: atom with more than 1 ring-stereo atoms found."
                    << std::endl;
              }
            }
          }
        } else {
          v = atom->getChiralTag();
        }
        (*atomCodes)[i] |= v << 18;  // 2 bits
      }
      if (0 != (CF_ATOM_AROMATIC & flags))
        (*atomCodes)[i] |= (atom->getIsAromatic() ? 1 : 0) << 20;  // 1 bit
      // if(0!=( & flags))
      //  (*atomCodes)[i] |= (atom-()) << 21; // 3 bits reserved
      if (0 != (CF_ISOTOPE & flags))
        (*atomCodes)[i] |= (atom->getIsotope()) << 24;
    }
  }

  if (bondCodes) {
    std::map<unsigned, bool> bondsInRing;
    const RingInfo::VECT_INT_VECT &rings = mol.getRingInfo()->bondRings();
    for (const auto &ring : rings)
      for (int b : ring)
        if (bondsInRing.end() == bondsInRing.find(b))
          bondsInRing[(unsigned)b] = true;

    unsigned n = mol.getNumBonds();
    bondCodes->resize(n);
    for (unsigned i = 0; i < n; i++) {
      if (0 == (CF_BOND_ALL & flags))  // NO LABELS
      {
        (*bondCodes)[i] = 1;
        continue;
      }
      const Bond *bond = mol.getBondWithIdx(i);
      (*bondCodes)[i] = 0;
      if (0 != (CF_BOND_ORDER & flags)) {
        unsigned order = bond->getBondType();
        if (0 == (CF_BOND_AROMATIZATION & flags))  // ignore aromatization
        {
          static const unsigned orderMatch[Bond::ZERO + 1] = {
              Bond::UNSPECIFIED, Bond::SINGLE,    Bond::DOUBLE,
              Bond::TRIPLE,      Bond::QUADRUPLE, Bond::QUINTUPLE,
              Bond::HEXTUPLE,
              Bond::SINGLE,     // ONEANDAHALF,
              Bond::DOUBLE,     // TWOANDAHALF,
              Bond::TRIPLE,     // THREEANDAHALF,
              Bond::QUADRUPLE,  // FOURANDAHALF,
              Bond::QUINTUPLE,  // FIVEANDAHALF,
              Bond::SINGLE,     // AROMATIC,
              Bond::IONIC,       Bond::HYDROGEN,  Bond::THREECENTER,
              Bond::DATIVEONE,   Bond::DATIVE,    Bond::DATIVEL,
              Bond::DATIVER,     Bond::OTHER,     Bond::ZERO};
          order = orderMatch[order];
        }
        (*bondCodes)[i] |= order;
      }
      if (0 != (CF_BOND_AROMATIZATION & flags))
        (*bondCodes)[i] |= ((bond->getIsAromatic() ? 1 : 0)) << 8;
      if (0 != (CF_BOND_CHIRALITY & flags)) {
        (*bondCodes)[i] |= bond->getStereo() << 9;
      }
      if (0 != (CF_BOND_IN_RING & flags))
        (*bondCodes)[i] |=
            (bondsInRing.end() != bondsInRing.find(bond->getIdx()) ? 1 : 0)
            << 11;  // 1 bit
    }
  }
}

//=============================================================================

HashCodeType generateMoleculeHashCode(
    const ROMol &mol, const std::vector<unsigned> *atomsToUse,
    const std::vector<unsigned> *bondsToUse,
    const std::vector<std::uint32_t> *atomCodes,
    const std::vector<std::uint32_t> *bondCodes) {
  MolFragment m;
  prepareMolFragment(m, mol, atomsToUse, bondsToUse);
  if (0 == m.getNumAtoms() || 0 == m.getNumBonds()) return 0;
  std::vector<std::uint32_t> atomLabels;
  std::vector<std::uint32_t> bondLabels;
  prepareLabels(atomLabels, bondLabels, mol, m, atomCodes, bondCodes);
  return computeMorganCodeHash(m, atomLabels, bondLabels);
}

void generateMoleculeHashSet(const ROMol &mol, HashSet &res,
                             const std::vector<unsigned> *atomsToUse,
                             const std::vector<unsigned> *bondsToUse) {
  res.Version = 100;  // v. 1.0
  res.Reserved = 0;

  MolFragment m;
  prepareMolFragment(m, mol, atomsToUse, bondsToUse);

  res.NumAtoms = m.getNumAtoms();
  res.NumBonds = m.getNumBonds();
  if (0 == m.getNumAtoms() || 0 == m.getNumBonds()) return;

  std::string formula = RDKit::Descriptors::calcMolFormula(mol);
  res.FormulaCRC32 = computeCRC32(formula.c_str(), formula.length());

  boost::uint64_t flags = 0;  // CodeFlags constants combination
  std::vector<std::uint32_t> atomCodes;
  std::vector<std::uint32_t> bondCodes;
  std::vector<std::uint32_t> atomLabels;
  std::vector<std::uint32_t> bondLabels;

  //        flags = CF_ATOM_ALL &(~(CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY |
  //        CF_ISOTOPE));
  flags = CF_ELEMENT | CF_CHARGE | CF_ATOM_AROMATIC;  /// | CF_VALENCE
  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);
  prepareLabels(atomLabels, bondLabels, mol, m, &atomCodes, &bondCodes);
  res.NonChiralAtomsHash = computeMorganCodeHash(m, atomLabels, bondLabels);

  flags = CF_BOND_ALL & (~CF_BOND_CHIRALITY);
  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);
  prepareLabels(atomLabels, bondLabels, mol, m, &atomCodes, &bondCodes);
  res.NonChiralBondsHash = computeMorganCodeHash(m, atomLabels, bondLabels);

  flags = CF_ATOM_CHIRALITY | CF_ISOTOPE;
  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);
  prepareLabels(atomLabels, bondLabels, mol, m, &atomCodes, &bondCodes);
  res.ChiralAtomsHash = computeMorganCodeHash(m, atomLabels, bondLabels);

  flags = CF_BOND_CHIRALITY;
  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);
  prepareLabels(atomLabels, bondLabels, mol, m, &atomCodes, &bondCodes);
  res.ChiralBondsHash = computeMorganCodeHash(m, atomLabels, bondLabels);

  flags = CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE;
  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);
  prepareLabels(atomLabels, bondLabels, mol, m, &atomCodes, &bondCodes);
  res.ChiralityHash = computeMorganCodeHash(m, atomLabels, bondLabels);
}

//=============================================================================
std::string generateMoleculeHashSet(const ROMol &mol,
                                    const std::vector<unsigned> *atomsToUse,
                                    const std::vector<unsigned> *bondsToUse) {
  std::string str;
  HashSet res;
  generateMoleculeHashSet(mol, res, atomsToUse, bondsToUse);
  // char buf[64];
  // snprintf(buf, sizeof(buf),"%u-%u-%u-", res.Version,
  // res.NumAtoms,res.NumBonds);
  str = (boost::format("%u-%u-%u-") % (res.Version) % (res.NumAtoms) %
         (res.NumBonds)).str();
  str += encode(&res.FormulaCRC32, sizeof(res.FormulaCRC32));
  str += "-";
  str += encode(&res.NonChiralAtomsHash, sizeof(res.NonChiralAtomsHash));
  str += "-";
  str += encode(&res.NonChiralBondsHash, sizeof(res.NonChiralBondsHash));
  str += "-";
  str += encode(&res.ChiralAtomsHash, sizeof(res.ChiralAtomsHash));
  str += "-";
  str += encode(&res.ChiralBondsHash, sizeof(res.ChiralBondsHash));
  str += "-";
  str += encode(&res.ChiralityHash, sizeof(res.ChiralityHash));
  return str.c_str();
}
//=============================================================================
// INTERNAL FUNCTIONS:
//=============================================================================
static HashCodeType computeMorganCodeHash(
    const MolFragment &mol, const std::vector<std::uint32_t> &atomLabels,
    const std::vector<std::uint32_t> &bondLabels) {
  size_t nv = mol.getNumAtoms();
  size_t ne = mol.getNumBonds();
  std::vector<HashCodeType> currCodes(nv);
  std::vector<HashCodeType> prevCodes(nv);
  size_t nIterations = mol.getNumBonds();
  if (nIterations > 5) nIterations = 5;

  for (unsigned molAtomIdx = 0; molAtomIdx < mol.getNumAtoms(); molAtomIdx++)
    currCodes[molAtomIdx] = atomLabels[mol.AtomsIdx[molAtomIdx]];

  for (size_t iter = 0; iter < nIterations; iter++) {
    for (size_t i = 0; i < nv; i++) prevCodes[i] = currCodes[i];

    for (size_t molBondIdx = 0; molBondIdx < ne; molBondIdx++) {
      const Bond *bond = mol.Bonds[molBondIdx];
      unsigned order = bondLabels[mol.BondsIdx[molBondIdx]];
      unsigned atom1 = mol.MolAtomIdxMap.find(bond->getBeginAtomIdx())->second;
      unsigned atom2 = mol.MolAtomIdxMap.find(bond->getEndAtomIdx())->second;
      std::uint32_t v1 = prevCodes[atom1];
      std::uint32_t v2 = prevCodes[atom2];

      currCodes[atom1] += v2 * v2 + (v2 + 23) * (order + 1721);
      currCodes[atom2] += v1 * v1 + (v1 + 23) * (order + 1721);
    }
  }

  HashCodeType result = 0;
  for (unsigned molAtomIdx = 0; molAtomIdx < nv; molAtomIdx++) {
    HashCodeType code = currCodes[molAtomIdx];
    result += code * (code + 6849) + 29;
  }
  return result;
}
//=============================================================================
static void prepareMolFragment(MolFragment &m, const ROMol &mol,
                               const std::vector<unsigned> *atomsToUse,
                               const std::vector<unsigned> *bondsToUse) {
  if (nullptr != atomsToUse && atomsToUse->empty()) atomsToUse = nullptr;
  if (nullptr != bondsToUse && bondsToUse->empty()) bondsToUse = nullptr;

  if (nullptr == atomsToUse && nullptr == bondsToUse)  // whole molecule
  {
    unsigned n = mol.getNumAtoms();
    m.AtomsIdx.resize(n);
    for (unsigned i = 0; i < n; i++) m.AtomsIdx[i] = i;

    n = mol.getNumBonds();
    m.BondsIdx.resize(n);
    for (unsigned i = 0; i < n; i++) m.BondsIdx[i] = i;
  } else if (nullptr !=
             atomsToUse)  // selected atoms only and all/selected bonds
                          // between them
  {
    std::map<unsigned, unsigned> addedBonds;
    unsigned n = atomsToUse->size();
    m.AtomsIdx.resize(n);
    for (unsigned i = 0; i < n; i++)  // add all selected atoms at first
      m.AtomsIdx[i] = (*atomsToUse)[i];

    for (unsigned i = 0; i < n; i++)  // add bonds between all selected atoms
    {
      ROMol::OEDGE_ITER beg, end;
      for (boost::tie(beg, end) =
               mol.getAtomBonds(mol.getAtomWithIdx(m.AtomsIdx[i]));
           beg != end; beg++) {
        const Bond *bond = &*((mol)[*beg]);
        if (addedBonds.end() != addedBonds.find(bond->getIdx()))
          continue;  // the bond has been already added
        if (nullptr != bondsToUse &&
            bondsToUse->end() ==
                find(bondsToUse->begin(), bondsToUse->end(), bond->getIdx()))
          continue;  // skip unselected bond

        unsigned endAtoms[2];
        endAtoms[0] = bond->getBeginAtomIdx();
        endAtoms[1] = bond->getEndAtomIdx();
        for (unsigned ai = 0; ai < 2 && nullptr != bond;
             ai++)  // both ending bonds of the atom
        {
          if (nullptr != atomsToUse &&
              atomsToUse->end() ==
                  find(atomsToUse->begin(), atomsToUse->end(), endAtoms[ai]))
            bond = nullptr;  // check if both ending atoms of the bond are
                             // selected by
                             // atoms filter
        }
        if (nullptr != bond) {
          addedBonds[bond->getIdx()] = m.BondsIdx.size();
          m.BondsIdx.push_back(bond->getIdx());
        }
      }
    }
  } else if (nullptr != bondsToUse)  // note that 0==atomsToUse in this case
  {
    std::map<unsigned, unsigned> addedAtoms;
    unsigned n = bondsToUse->size();
    m.BondsIdx.resize(n);
    for (unsigned i = 0; i < n; i++) {
      const Bond *bond = mol.getBondWithIdx(i);
      m.BondsIdx[i] = bond->getIdx();

      unsigned endAtoms[2];
      endAtoms[0] = bond->getBeginAtomIdx();
      endAtoms[1] = bond->getEndAtomIdx();
      for (unsigned ai = 0; ai < 2 && nullptr != bond;
           ai++)  // both ending bonds of the atom
      {
        if (addedAtoms.end() == addedAtoms.find(endAtoms[ai])) {
          addedAtoms[endAtoms[ai]] = addedAtoms.size();
          m.AtomsIdx.push_back(
              endAtoms[ai]);  // the atom has NOT been already added
        }
      }
    }
  }

  unsigned n;
  n = m.getNumAtoms();
  m.Atoms.resize(n);
  for (unsigned i = 0; i < n; i++) {
    m.Atoms[i] = mol.getAtomWithIdx(m.AtomsIdx[i]);
    m.MolAtomIdxMap[m.AtomsIdx[i]] = i;
  }

  n = m.getNumBonds();
  m.Bonds.resize(n);
  for (unsigned i = 0; i < n; i++) {
    m.Bonds[i] = mol.getBondWithIdx(m.BondsIdx[i]);
  }
}
//=============================================================================
static void prepareLabels(std::vector<std::uint32_t> &atomLabels,
                          std::vector<std::uint32_t> &bondLabels,
                          const ROMol &mol, const MolFragment &m,
                          const std::vector<std::uint32_t> *atomCodes,
                          const std::vector<std::uint32_t> *bondCodes) {
  RDUNUSED_PARAM(mol);
  unsigned n;
  n = m.getNumAtoms();
  atomLabels.resize(n);
  for (unsigned i = 0; i < n; i++) {
    atomLabels[i] = atomCodes ? (*atomCodes)[m.AtomsIdx[i]] : 1;
  }

  n = m.getNumBonds();
  bondLabels.resize(n);
  for (unsigned i = 0; i < n; i++) {
    bondLabels[i] = bondCodes ? (*bondCodes)[m.BondsIdx[i]] : 1;
  }
}
//=============================================================================
}
}
