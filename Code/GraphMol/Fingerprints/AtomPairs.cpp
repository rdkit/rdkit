//
//  Copyright (C) 2007-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>

namespace RDKit {
namespace AtomPairs {

template <typename T1, typename T2>
void updateElement(SparseIntVect<T1> &v, T2 elem) {
  v.setVal(elem, v.getVal(elem) + 1);
}

template <typename T1>
void updateElement(ExplicitBitVect &v, T1 elem) {
  v.setBit(elem % v.getNumBits());
}

template <typename T>
void setAtomPairBit(std::uint32_t i, std::uint32_t j, std::uint32_t nAtoms,
                    const std::vector<std::uint32_t> &atomCodes,
                    const double *dm, T *bv, unsigned int minLength,
                    unsigned int maxLength, bool includeChirality) {
  auto dist = static_cast<unsigned int>(floor(dm[i * nAtoms + j]));
  if (dist >= minLength && dist <= maxLength) {
    std::uint32_t bitId =
        getAtomPairCode(atomCodes[i], atomCodes[j], dist, includeChirality);
    updateElement(*bv, static_cast<std::uint32_t>(bitId));
  }
}

SparseIntVect<std::int32_t> *getAtomPairFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality,
    bool use2D, int confId) {
  return getAtomPairFingerprint(mol, 1, maxPathLen - 1, fromAtoms, ignoreAtoms,
                                atomInvariants, includeChirality, use2D,
                                confId);
};

SparseIntVect<std::int32_t> *getAtomPairFingerprint(
    const ROMol &mol, unsigned int minLength, unsigned int maxLength,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality,
    bool use2D, int confId) {
  PRECONDITION(minLength <= maxLength, "bad lengths provided");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }

  auto *res = new SparseIntVect<std::int32_t>(
      1 << (numAtomPairFingerprintBits + 2 * (includeChirality ? 2 : 0)));
  const double *dm;
  if (use2D) {
    dm = MolOps::getDistanceMat(*lmol);
  } else {
    dm = MolOps::get3DDistanceMat(*lmol, confId);
  }
  const unsigned int nAtoms = lmol->getNumAtoms();

  std::vector<std::uint32_t> atomCodes;
  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    } else {
      atomCodes.push_back((*atomInvariants)[(*atomItI)->getIdx()] %
                          ((1 << codeSize) - 1));
    }
  }

  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    unsigned int i = (*atomItI)->getIdx();
    if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
                           ignoreAtoms->end()) {
      continue;
    }
    if (!fromAtoms) {
      for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
           atomItJ != lmol->endAtoms(); ++atomItJ) {
        unsigned int j = (*atomItJ)->getIdx();
        if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                     j) != ignoreAtoms->end()) {
          continue;
        }
        setAtomPairBit(i, j, nAtoms, atomCodes, dm, res, minLength, maxLength,
                       includeChirality);
      }
    } else {
      for (auto j : *fromAtoms) {
        if (j != i) {
          if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                       j) != ignoreAtoms->end()) {
            continue;
          }
          setAtomPairBit(i, j, nAtoms, atomCodes, dm, res, minLength, maxLength,
                         includeChirality);
        }
      }
    }
  }
  return res;
}

SparseIntVect<std::int32_t> *getHashedAtomPairFingerprint(
    const ROMol &mol, unsigned int nBits, unsigned int minLength,
    unsigned int maxLength, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality,
    bool use2D, int confId) {
  PRECONDITION(minLength <= maxLength, "bad lengths provided");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }

  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
      RDKit::AtomPair::getAtomPairGenerator<std::uint32_t>(
          minLength, maxLength, includeChirality, use2D, nullptr, true, nBits)};
  FingerprintFuncArguments args;
  args.fromAtoms = fromAtoms;
  args.ignoreAtoms = ignoreAtoms;
  args.customAtomInvariants = atomInvariants;
  args.confId = confId;
  auto siv = fpgen->getCountFingerprint(*lmol, args);
  auto *res = new SparseIntVect<std::int32_t>(nBits);
  for (auto v : siv->getNonzeroElements()) {
    res->setVal(v.first, v.second);
  }
  return res;
}

ExplicitBitVect *getHashedAtomPairFingerprintAsBitVect(
    const ROMol &mol, unsigned int nBits, unsigned int minLength,
    unsigned int maxLength, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants,
    unsigned int nBitsPerEntry, bool includeChirality, bool use2D, int confId) {
  PRECONDITION(minLength <= maxLength, "bad lengths provided");
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  static int bounds[4] = {1, 2, 4, 8};

  unsigned int blockLength = nBits / nBitsPerEntry;
  SparseIntVect<std::int32_t> *sres = getHashedAtomPairFingerprint(
      mol, blockLength, minLength, maxLength, fromAtoms, ignoreAtoms,
      atomInvariants, includeChirality, use2D, confId);
  auto *res = new ExplicitBitVect(nBits);
  if (nBitsPerEntry != 4) {
    for (auto val : sres->getNonzeroElements()) {
      for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
        if (val.second > static_cast<int>(i)) {
          res->setBit(val.first * nBitsPerEntry + i);
        }
      }
    }
  } else {
    for (auto val : sres->getNonzeroElements()) {
      for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
        if (val.second >= bounds[i]) {
          res->setBit(val.first * nBitsPerEntry + i);
        }
      }
    }
  }
  delete sres;
  return res;
}

SparseIntVect<boost::int64_t> *getTopologicalTorsionFingerprint(
    const ROMol &mol, unsigned int targetSize,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }
  boost::uint64_t sz = 1;
  sz = (sz << (targetSize *
               (codeSize + (includeChirality ? numChiralBits : 0))));
  // NOTE: this -1 is incorrect but it's needed for backwards compatibility.
  //  hopefully we'll never have a case with a torsion that hits this.
  //
  //  mmm, bug compatible.
  sz -= 1;
  auto *res = new SparseIntVect<boost::int64_t>(sz);

  std::vector<std::uint32_t> atomCodes;
  atomCodes.reserve(lmol->getNumAtoms());
  for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
       atomItI != lmol->endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    } else {
      // need to add to the atomCode here because we subtract off up to 2 below
      // as part of the branch correction
      atomCodes.push_back(
          (*atomInvariants)[(*atomItI)->getIdx()] % ((1 << codeSize) - 1) + 2);
    }
  }

  boost::dynamic_bitset<> *fromAtomsBV = nullptr;
  if (fromAtoms) {
    fromAtomsBV = new boost::dynamic_bitset<>(lmol->getNumAtoms());
    for (auto fAt : *fromAtoms) {
      fromAtomsBV->set(fAt);
    }
  }
  boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
  if (ignoreAtoms) {
    ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    for (auto fAt : *ignoreAtoms) {
      ignoreAtomsBV->set(fAt);
    }
  }
  boost::dynamic_bitset<> pAtoms(lmol->getNumAtoms());
  PATH_LIST paths = findAllPathsOfLengthN(*lmol, targetSize, false);
  for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
       ++pathIt) {
    bool keepIt = true;
    if (fromAtomsBV) {
      keepIt = false;
    }
    std::vector<std::uint32_t> pathCodes;
    const PATH_TYPE &path = *pathIt;
    if (fromAtomsBV) {
      if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
          fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
        keepIt = true;
      }
    }
    if (keepIt && ignoreAtomsBV) {
      for (auto pElem : path) {
        if (ignoreAtomsBV->test(pElem)) {
          keepIt = false;
          break;
        }
      }
    }
    if (keepIt) {
      pAtoms.reset();
      for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
        // look for a cycle that doesn't start at the first atom
        // we can't effectively canonicalize these at the moment
        // (was github #811)
        if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
          pathCodes.clear();
          break;
        }
        pAtoms.set(*pIt);
        unsigned int code = atomCodes[*pIt] - 1;
        // subtract off the branching number:
        if (pIt != path.begin() && pIt + 1 != path.end()) {
          --code;
        }
        pathCodes.push_back(code);
      }
      if (pathCodes.size()) {
        boost::int64_t code =
            getTopologicalTorsionCode(pathCodes, includeChirality);
        updateElement(*res, code);
      }
    }
  }
  delete fromAtomsBV;
  delete ignoreAtomsBV;

  return res;
}

namespace {
template <typename T>
void TorsionFpCalc(T *res, const ROMol &mol, unsigned int nBits,
                   unsigned int targetSize,
                   const std::vector<std::uint32_t> *fromAtoms,
                   const std::vector<std::uint32_t> *ignoreAtoms,
                   const std::vector<std::uint32_t> *atomInvariants,
                   bool includeChirality) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen{
      RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
          includeChirality, targetSize, nullptr, true, nBits)};
  FingerprintFuncArguments args;
  args.fromAtoms = fromAtoms;
  args.ignoreAtoms = ignoreAtoms;
  args.customAtomInvariants = atomInvariants;
  auto siv = fpgen->getCountFingerprint(*lmol, args);
  for (auto v : siv->getNonzeroElements()) {
    res->setVal(v.first, v.second);
  }
}
}  // namespace
SparseIntVect<boost::int64_t> *getHashedTopologicalTorsionFingerprint(
    const ROMol &mol, unsigned int nBits, unsigned int targetSize,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants, bool includeChirality) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  auto *res = new SparseIntVect<boost::int64_t>(nBits);
  TorsionFpCalc(res, mol, nBits, targetSize, fromAtoms, ignoreAtoms,
                atomInvariants, includeChirality);
  return res;
}

ExplicitBitVect *getHashedTopologicalTorsionFingerprintAsBitVect(
    const ROMol &mol, unsigned int nBits, unsigned int targetSize,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const std::vector<std::uint32_t> *atomInvariants,
    unsigned int nBitsPerEntry, bool includeChirality) {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  static int bounds[4] = {1, 2, 4, 8};
  unsigned int blockLength = nBits / nBitsPerEntry;
  auto *sres = new SparseIntVect<boost::int64_t>(blockLength);
  TorsionFpCalc(sres, mol, blockLength, targetSize, fromAtoms, ignoreAtoms,
                atomInvariants, includeChirality);
  auto *res = new ExplicitBitVect(nBits);

  if (nBitsPerEntry != 4) {
    for (auto val : sres->getNonzeroElements()) {
      for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
        if (val.second > static_cast<int>(i)) {
          res->setBit(val.first * nBitsPerEntry + i);
        }
      }
    }
  } else {
    for (auto val : sres->getNonzeroElements()) {
      for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
        if (val.second >= bounds[i]) {
          res->setBit(val.first * nBitsPerEntry + i);
        }
      }
    }
  }
  delete sres;
  return res;
}
}  // end of namespace AtomPairs
}  // end of namespace RDKit
