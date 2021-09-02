//
//  Copyright (c) 2009-2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Created by Greg Landrum, July 2008
//
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <algorithm>

#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <RDGeneral/Exceptions.h>

namespace {
class ss_matcher {
 public:
  ss_matcher(){};
  ss_matcher(const std::string &pattern) {
    RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    TEST_ASSERT(p);
    m_matcher.reset(p);
  };

  // const RDKit::ROMOL_SPTR &getMatcher() const { return m_matcher; };
  const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };

 private:
  RDKit::ROMOL_SPTR m_matcher;
};
}  // namespace

namespace RDKit {
namespace MorganFingerprints {

uint32_t updateElement(SparseIntVect<uint32_t> &v, unsigned int elem,
                       bool counting) {
  uint32_t bit = elem % v.getLength();
  if (counting) {
    v.setVal(bit, v.getVal(bit) + 1);
  } else {
    v.setVal(bit, 1);
  }
  return bit;
}
uint32_t updateElement(ExplicitBitVect &v, unsigned int elem, bool) {
  uint32_t bit = elem % v.getNumBits();
  v.setBit(bit);
  return bit;
}

template <typename T>
void calcFingerprint(const ROMol &mol, unsigned int radius,
                     std::vector<uint32_t> *invariants,
                     const std::vector<uint32_t> *fromAtoms, bool useChirality,
                     bool useBondTypes, bool useCounts,
                     bool onlyNonzeroInvariants, BitInfoMap *atomsSettingBits,
                     bool includeRedundantEnvironments, T &res) {
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (useChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }
  unsigned int nAtoms = lmol->getNumAtoms();
  bool owner = false;
  if (!invariants) {
    invariants = new std::vector<uint32_t>(nAtoms);
    owner = true;
    getConnectivityInvariants(*lmol, *invariants);
  }
  // Make a copy of the invariants:
  std::vector<uint32_t> invariantCpy(nAtoms);
  std::copy(invariants->begin(), invariants->end(), invariantCpy.begin());

  // add the round 0 invariants to the result:
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (!fromAtoms || std::find(fromAtoms->begin(), fromAtoms->end(), i) !=
                          fromAtoms->end()) {
      if (!onlyNonzeroInvariants || (*invariants)[i]) {
        uint32_t bit = updateElement(res, (*invariants)[i], useCounts);
        if (atomsSettingBits) {
          (*atomsSettingBits)[bit].push_back(std::make_pair(i, 0));
        }
      }
    }
  }

  boost::dynamic_bitset<> chiralAtoms(nAtoms);

  // these are the neighborhoods that have already been added
  // to the fingerprint
  std::vector<boost::dynamic_bitset<>> neighborhoods;
  // these are the environments around each atom:
  std::vector<boost::dynamic_bitset<>> atomNeighborhoods(
      nAtoms, boost::dynamic_bitset<>(mol.getNumBonds()));
  boost::dynamic_bitset<> deadAtoms(nAtoms);

  boost::dynamic_bitset<> includeAtoms(nAtoms);
  if (fromAtoms) {
    for (auto idx : *fromAtoms) {
      includeAtoms.set(idx, 1);
    }
  } else {
    includeAtoms.set();
  }

  std::vector<unsigned int> atomOrder(nAtoms);
  if (onlyNonzeroInvariants) {
    std::vector<std::pair<int32_t, uint32_t>> ordering;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (!(*invariants)[i]) {
        ordering.emplace_back(1, i);
      } else {
        ordering.emplace_back(0, i);
      }
    }
    std::sort(ordering.begin(), ordering.end());
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = ordering[i].second;
    }
  } else {
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = i;
    }
  }
  // now do our subsequent rounds:
  for (unsigned int layer = 0; layer < radius; ++layer) {
    std::vector<uint32_t> roundInvariants(nAtoms);
    std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
        atomNeighborhoods;
    std::vector<AccumTuple> neighborhoodsThisRound;

    for (auto atomIdx : atomOrder) {
      if (!deadAtoms[atomIdx]) {
        const Atom *tAtom = lmol->getAtomWithIdx(atomIdx);
        if (!tAtom->getDegree()) {
          deadAtoms.set(atomIdx, 1);
          continue;
        }
        std::vector<std::pair<int32_t, uint32_t>> nbrs;
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = lmol->getAtomBonds(tAtom);
        while (beg != end) {
          const Bond *bond = mol[*beg];
          roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;

          unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
          roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

          int32_t bt = 1;
          if (useBondTypes) {
            if (!useChirality || bond->getBondType() != Bond::DOUBLE ||
                bond->getStereo() == Bond::STEREONONE) {
              bt = static_cast<int32_t>(bond->getBondType());
            } else {
              const int32_t stereoOffset = 100;
              const int32_t bondTypeOffset = 10;
              bt = stereoOffset +
                   bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
                   static_cast<int32_t>(bond->getStereo());
            }
          }
          nbrs.emplace_back(bt, (*invariants)[oIdx]);

          ++beg;
        }

        // sort the neighbor list:
        std::sort(nbrs.begin(), nbrs.end());
        // and now calculate the new invariant and test if the atom is newly
        // "chiral"
        std::uint32_t invar = layer;
        gboost::hash_combine(invar, (*invariants)[atomIdx]);
        bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
        for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
                 nbrs.begin();
             it != nbrs.end(); ++it) {
          // add the contribution to the new invariant:
          gboost::hash_combine(invar, *it);

          // std::cerr<<"     "<<atomIdx<<": "<<it->first<<" "<<it->second<<" ->
          // "<<invar<<std::endl;

          // update our "chirality":
          if (useChirality && looksChiral && chiralAtoms[atomIdx]) {
            if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
              looksChiral = false;
            } else if (it != nbrs.begin() && it->second == (it - 1)->second) {
              looksChiral = false;
            }
          }
        }
        if (useChirality && looksChiral) {
          chiralAtoms[atomIdx] = 1;
          // add an extra value to the invariant to reflect chirality:
          std::string cip = "";
          tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
          if (cip == "R") {
            gboost::hash_combine(invar, 3);
          } else if (cip == "S") {
            gboost::hash_combine(invar, 2);
          } else {
            gboost::hash_combine(invar, 1);
          }
        }
        roundInvariants[atomIdx] = static_cast<uint32_t>(invar);
        neighborhoodsThisRound.push_back(
            boost::make_tuple(roundAtomNeighborhoods[atomIdx],
                              static_cast<uint32_t>(invar), atomIdx));
        if (!includeRedundantEnvironments &&
            std::find(neighborhoods.begin(), neighborhoods.end(),
                      roundAtomNeighborhoods[atomIdx]) != neighborhoods.end()) {
          // we have seen this exact environment before, this atom
          // is now out of consideration:
          deadAtoms[atomIdx] = 1;
          // std::cerr<<"   atom: "<< atomIdx <<" is dead."<<std::endl;
        }
      }
    }
    std::sort(neighborhoodsThisRound.begin(), neighborhoodsThisRound.end());
    for (std::vector<AccumTuple>::const_iterator iter =
             neighborhoodsThisRound.begin();
         iter != neighborhoodsThisRound.end(); ++iter) {
      // if we haven't seen this exact environment before, update the
      // fingerprint:
      if (includeRedundantEnvironments ||
          std::find(neighborhoods.begin(), neighborhoods.end(),
                    iter->get<0>()) == neighborhoods.end()) {
        if (!onlyNonzeroInvariants || invariantCpy[iter->get<2>()]) {
          if (includeAtoms[iter->get<2>()]) {
            uint32_t bit = updateElement(res, iter->get<1>(), useCounts);
            if (atomsSettingBits) {
              (*atomsSettingBits)[bit].push_back(
                  std::make_pair(iter->get<2>(), layer + 1));
            }
          }
          if (!fromAtoms || std::find(fromAtoms->begin(), fromAtoms->end(),
                                      iter->get<2>()) != fromAtoms->end()) {
            neighborhoods.push_back(iter->get<0>());
          }
        }
        // std::cerr<<" layer: "<<layer<<" atom: "<<iter->get<2>()<<" "
        // <<iter->get<0>()<< " " << iter->get<1>() << " " <<
        // deadAtoms[iter->get<2>()]<<std::endl;
      } else {
        // we have seen this exact environment before, this atom
        // is now out of consideration:
        // std::cerr<<"   atom: "<< iter->get<2>()<<" is dead."<<std::endl;
        deadAtoms[iter->get<2>()] = 1;
      }
    }

    // the invariants from this round become the global invariants:
    std::copy(roundInvariants.begin(), roundInvariants.end(),
              invariants->begin());

    atomNeighborhoods = roundAtomNeighborhoods;
  }

  if (owner) {
    delete invariants;
  }
}

SparseIntVect<uint32_t> *getFingerprint(
    const ROMol &mol, unsigned int radius, std::vector<uint32_t> *invariants,
    const std::vector<uint32_t> *fromAtoms, bool useChirality,
    bool useBondTypes, bool useCounts, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  SparseIntVect<uint32_t> *res;
  res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());
  calcFingerprint(mol, radius, invariants, fromAtoms, useChirality,
                  useBondTypes, useCounts, onlyNonzeroInvariants,
                  atomsSettingBits, includeRedundantEnvironments, *res);
  return res;
}
SparseIntVect<uint32_t> *getHashedFingerprint(
    const ROMol &mol, unsigned int radius, unsigned int nBits,
    std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  if (nBits == 0) {
    throw ValueErrorException("nBits can not be zero");
  }
  SparseIntVect<uint32_t> *res;
  res = new SparseIntVect<uint32_t>(nBits);
  calcFingerprint(mol, radius, invariants, fromAtoms, useChirality,
                  useBondTypes, true, onlyNonzeroInvariants, atomsSettingBits,
                  includeRedundantEnvironments, *res);
  return res;
}

ExplicitBitVect *getFingerprintAsBitVect(
    const ROMol &mol, unsigned int radius, unsigned int nBits,
    std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  auto *res = new ExplicitBitVect(nBits);
  calcFingerprint(mol, radius, invariants, fromAtoms, useChirality,
                  useBondTypes, false, onlyNonzeroInvariants, atomsSettingBits,
                  includeRedundantEnvironments, *res);
  return res;
}

}  // end of namespace MorganFingerprints
}  // end of namespace RDKit
