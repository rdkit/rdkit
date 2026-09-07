//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <DistGeom/ZMatrix.h>
#include "ZMatrixBuilder.h"
#include "BoundsMatrixBuilderDetails.h"
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <ForceField/UFF/BondStretch.h>
#include <Geometry/Utils.h>

#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <Numerics/SymmMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/multiprecision/integer.hpp>
#include <optional>
#include <ranges>
#include <DistGeom/ChiralSet.h>
#include <vector>

namespace RDKit {
namespace DGeomHelpers {

using Type14References = std::vector<std::optional<std::tuple<
    std::optional<unsigned int>, unsigned int, std::optional<unsigned int>>>>;

struct StackElem {
  unsigned int atomIdx, precursorIdx;
};

void addNeighborsToStack(const unsigned int atomIdx,
                         const unsigned int precursorIdx, const ROMol &mol,
                         std::vector<StackElem> &stack,
                         const boost::dynamic_bitset<> &visitedAtoms,
                         const DistGeom::ZMatrix &zmat,
                         const InternalCoordinates &internalCoords) {
  std::vector<const RDKit::Bond *> bonds;
  std::vector<unsigned int> ringClosures;

  for (const auto &bnd : mol.atomBonds(mol.getAtomWithIdx(atomIdx))) {
    const auto nbratomIdx = bnd->getOtherAtomIdx(atomIdx);
    if (nbratomIdx == precursorIdx) {
      continue;
    }
    if (visitedAtoms[nbratomIdx]) {
      ringClosures.push_back(nbratomIdx);
    } else {
      bonds.push_back(bnd);
    }
  }

  if (bonds.empty() && ringClosures.empty()) {
    return;
  }

  const auto &references = zmat.getReferences(atomIdx);

  if (references.angleRef && !bonds.empty()) {
    // if torsions are given, we want to add the most constraint one at last to
    // the stack
    const auto bondIdx1 =
        mol.getBondBetweenAtoms(*references.angleRef, precursorIdx)->getIdx();

    const auto bondIdx2 =
        mol.getBondBetweenAtoms(precursorIdx, atomIdx)->getIdx();

    const auto bondToTorsion =
        [&](const Bond *bond) -> const DistGeom::TorsionCandidates & {
      return internalCoords.torsionRange.at(
          getUnifiedId(bondIdx1, bondIdx2, bond->getIdx(), mol.getNumBonds()));
    };

    const auto first =
        std::ranges::min_element(bonds, DistGeom::less, bondToTorsion);

    const auto next = std::next(
        first);  // it should be handelt first => adding it as last to the stack
    std::ranges::rotate(bonds, next == bonds.end() ? bonds.begin() : next);
  }

  for (const auto *bond : bonds) {
    stack.emplace_back(bond->getOtherAtomIdx(atomIdx), atomIdx);
  }

  // ring closures should be handles first
  for (const auto nbrIdx : ringClosures) {
    stack.emplace_back(nbrIdx, atomIdx);
  }
}

void ringClosure(const unsigned int atomIdx, const unsigned int precursorIdx,
                 Type14References &references) {
  const auto &[tRef, ref1, ref2Opt] =
      references[precursorIdx]
          .value();  // there mus be a value, since precursor was already
                     // set (unless for the first atom which cannot be a ring
                     // closure!)

  std::optional<DistGeom::ZMatrix::TorsionDependence> torsiondependence;
  if (tRef) {
    // we only need this for the chiral correction but the offset
    torsiondependence = {0.0, tRef.value()};
  }

  if (ref2Opt && atomIdx == ref2Opt.value()) {
    // this is the case for small rins (3 atoms) or highly fused sysems
    return;
  }

  references[precursorIdx] = {{atomIdx}, ref1, ref2Opt};

  CHECK_INVARIANT(references[atomIdx],
                  "Ring Closure -> references must have been set before");
  {
    const auto &[_torsionDependence, _angleRef, _torsionRef] =
        *references[atomIdx];
    // note: this only makes a difference if the ring closure is more or
    // less
    // chemically valid -> we will still do it here for the case that we
    // have
    // perfect rings :)
    references[atomIdx] = {precursorIdx, _angleRef, _torsionRef};
  }
}

void addElement(const Atom *atom, const unsigned int precursorIdx,
                const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
                const InternalCoordinates &internalCoords,
                Type14References &references) {
  const auto &[tRef, ref1, ref2Opt] =
      references[precursorIdx]
          .value();  // there mus be a value since the precursor was already
  // set

  const unsigned int bndIdx1 =
      mol.getBondBetweenAtoms(atom->getIdx(), precursorIdx)->getIdx();

  const unsigned int bndIdx2 =
      mol.getBondBetweenAtoms(ref1, precursorIdx)->getIdx();

  std::optional<double> bl = internalCoords.lengths[bndIdx1];

  if (internalCoords.angles.find(static_cast<std::uint64_t>(
          getUnifiedId(bndIdx1, bndIdx2, mol.getNumBonds()))) ==
      internalCoords.angles.end()) {
    std::cerr << "Angle " << bndIdx1 << " " << bndIdx2 << " not found\n";
  }

  double ba = internalCoords.angles.at(static_cast<std::uint64_t>(
      getUnifiedId(bndIdx1, bndIdx2, mol.getNumBonds())));

  if (RDKit::feq(ba, M_PI, 1e-8)) {
    // to avoid colinearity
    ba += 1e-6;
  }

  std::optional<DistGeom::TorsionCandidates> torsion = std::nullopt;
  std::optional<DistGeom::ZMatrix::TorsionDependence> torsiondependence;
  if (tRef) {
    // this also means that we do not need a torsion reference
    const double offset =
        2 * M_PI /
        static_cast<double>(
            std::min(mol.getAtomWithIdx(precursorIdx)->getDegree(), 5u) - 1);
    torsiondependence = {offset, tRef.value()};
  } else if (ref2Opt) {  // only if a torsion reference exists (this does not
                         // hold for the 3rd element in the matrix therefore the
                         // check)
    if (internalCoords.torsionRange.find(
            static_cast<std::uint64_t>(getUnifiedId(
                bndIdx1, bndIdx2,
                mol.getBondBetweenAtoms(ref1, ref2Opt.value())->getIdx(),
                mol.getNumBonds()))) == internalCoords.torsionRange.end()) {
      std::cerr << "Torsion " << bndIdx1 << " " << bndIdx2 << " "
                << mol.getBondBetweenAtoms(ref1, ref2Opt.value())->getIdx()
                << " not found\n";
    }
    torsion = internalCoords.torsionRange.at(static_cast<std::uint64_t>(
        getUnifiedId(bndIdx1, bndIdx2,
                     mol.getBondBetweenAtoms(ref1, ref2Opt.value())->getIdx(),
                     mol.getNumBonds())));
  }

  zmat->addElement(atom->getIdx(), bl, precursorIdx, {ba}, {ref1}, torsion,
                   torsiondependence ? std::nullopt : ref2Opt,
                   torsiondependence);

  std::get<0>(references[precursorIdx].value()) = {atom->getIdx()};

  assert(!references[atom->getIdx()]);  // invariant since we visit every

  references[atom->getIdx()] = {{}, precursorIdx, {ref1}};
}

void setMoleculeDFS(const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
                    const InternalCoordinates &internalCoords) {
  unsigned int idx1, idx2;

  if (mol.getNumAtoms() == 1) {
    setMoleculeDFS(mol, zmat, internalCoords, 0, 0, 0);
    return;
  }

  if (mol.getNumAtoms() == 1) {
    setMoleculeDFS(mol, zmat, internalCoords, 0, 1, 0);
  }

  if (internalCoords.torsionRange.empty()) {
    // since connected component -> there must be at least one angle
    const auto startAngleId = internalCoords.angles.begin()->first;
    const auto _temp = unifiedIdToBondIds<2>(startAngleId, mol.getNumBonds());
    idx1 = _temp[0];
    idx2 = _temp[1];
  } else {
    // starting with most constraint one
    const auto startTorsionId = internalCoords.torsionRange.begin()->first;
    const auto _temp = unifiedIdToBondIds<3>(startTorsionId, mol.getNumBonds());
    idx1 = _temp[0];
    idx2 = _temp[1];
  }

  const auto b1 = mol.getBondWithIdx(idx1);
  const auto b2 = mol.getBondWithIdx(idx2);

  auto startAtm = b1->getBeginAtomIdx();
  auto atom2Idx = b1->getEndAtomIdx();

  if (b2->getBeginAtomIdx() == startAtm || b2->getEndAtomIdx() == startAtm) {
    std::swap(startAtm, atom2Idx);
  }

  const auto atom3Idx = b2->getOtherAtomIdx(atom2Idx);

  setMoleculeDFS(mol, zmat, internalCoords, startAtm, atom2Idx, atom3Idx);
}

void setMoleculeDFS(const ROMol &mol, std::shared_ptr<DistGeom::ZMatrix> zmat,
                    const InternalCoordinates &internalCoords,
                    const unsigned int startAtomIdx,
                    const unsigned int atomIdx2, const unsigned int atomIdx3) {
  // init bookkeeping structure(s) -> torsionReferences
  Type14References references(
      mol.getNumAtoms());  // similar to zmatrix but gets updated + tracks
                           // backwards => important for first atoms

  boost::dynamic_bitset<> visitedAtoms(mol.getNumAtoms());
  boost::dynamic_bitset<> visitedBonds(mol.getNumBonds());

  std::vector<StackElem> stack;
  stack.reserve(mol.getNumAtoms() / 2);  // just an approximation

  zmat->addElement(startAtomIdx);

  switch (mol.getNumAtoms()) {
    case 2:
      zmat->addElement(
          atomIdx2,
          internalCoords.lengths[mol.getBondBetweenAtoms(startAtomIdx, atomIdx2)
                                     ->getIdx()],
          startAtomIdx);
      [[fallthrough]];
    case 1:
      return;
  }

  addNeighborsToStack(startAtomIdx, startAtomIdx, mol, stack, visitedAtoms,
                      *zmat, internalCoords);
  visitedAtoms.set(startAtomIdx);

  zmat->addElement(
      atomIdx2,
      internalCoords
          .lengths[mol.getBondBetweenAtoms(startAtomIdx, atomIdx2)->getIdx()],
      startAtomIdx);

  addNeighborsToStack(atomIdx2, startAtomIdx, mol, stack, visitedAtoms, *zmat,
                      internalCoords);

  visitedAtoms.set(atomIdx2);
  visitedBonds.set(mol.getBondBetweenAtoms(atomIdx2, startAtomIdx)->getIdx());

  references[startAtomIdx] = {std::nullopt, atomIdx2, atomIdx3};
  references[atomIdx2] = {std::nullopt, startAtomIdx, std::nullopt};

  // add third prio again to stack to ensure that we start with a torsion
  if (mol.getBondBetweenAtoms(atomIdx2, atomIdx3)) {
    stack.emplace_back(atomIdx3, atomIdx2);
  } else {
    stack.emplace_back(atomIdx3, startAtomIdx);
  }

  while (stack.size()) {
    const auto &[idx, precursor] = stack.back();
    stack.pop_back();

    unsigned int bndIdx = mol.getBondBetweenAtoms(idx, precursor)->getIdx();

    if (visitedBonds[bndIdx]) {
      continue;
    }

    if (!visitedAtoms[idx]) {
      addElement(mol.getAtomWithIdx(idx), precursor, mol, zmat, internalCoords,
                 references);
      visitedAtoms.set(idx);

      addNeighborsToStack(idx, precursor, mol, stack, visitedAtoms, *zmat,
                          internalCoords);
    } else {
      ringClosure(idx, precursor, references);
    }

    visitedBonds.set(bndIdx);
  }
}

void correctChiralCenters(const ROMol &mol,
                          std::shared_ptr<DistGeom::ZMatrix> zmat) {
  boost::dynamic_bitset<> visited{mol.getNumAtoms()};

  for (const auto &row :
       *zmat |
           std::views::drop(
               4)  // do not visit first 4 element -> the four atoms can only
                   // span a torsion -> we cannot have a torsion dependence here
           | std::views::reverse) {
    const auto centerIdx = *row.internal.bondRef;
    if (visited[centerIdx]) {
      // we already dealt with this atom
      continue;
    }

    const auto *center = mol.getAtomWithIdx(centerIdx);
    const auto chiralTag = center->getChiralTag();
    if ((chiralTag != Atom::CHI_TETRAHEDRAL_CW &&  // only consider tetrahereal
         chiralTag != Atom::CHI_TETRAHEDRAL_CCW) ||
        !row.torsionDependence) {  // TODO maybe removenot representable
      visited.set(centerIdx);
      continue;
    }

    const auto axisatomIdx = *row.internal.angleRef;

    // collect bonds in reverse order of setting them
    auto refIdx = row.atomIdx;
    std::vector<unsigned int> dependentIdxs;
    while (zmat->getTorsionReference(refIdx)) {
      dependentIdxs.emplace_back(refIdx);
      refIdx = zmat->getTorsionReference(refIdx)->reference;
    }

    const auto anchorIdx = refIdx;  // the first one that we have placed

    const auto *bnd1 = mol.getBondBetweenAtoms(centerIdx, axisatomIdx);
    const auto *bnd2 = mol.getBondBetweenAtoms(centerIdx, anchorIdx);

    if (!bnd1 || !bnd2 || bnd1 == bnd2) {
      // invalid center due to ring closure or for fused systems
      // we cannot correct chirality here
      visited.set(centerIdx);
      continue;
    }

    INT_LIST currentPertOrder;
    currentPertOrder.emplace_back(bnd1->getIdx());
    currentPertOrder.emplace_back(bnd2->getIdx());

    for (unsigned int &dependentIdx : std::views::reverse(dependentIdxs)) {
      const auto *bnd = mol.getBondBetweenAtoms(centerIdx, dependentIdx);
      if (!bnd ||
          std::ranges::find(currentPertOrder, bnd->getIdx()) !=
              currentPertOrder.end()) {  // this can happen in fused systems
        break;
      }

      currentPertOrder.emplace_back(bnd->getIdx());
    }

    if (currentPertOrder.size() != center->getDegree()) {
      // this can happen for ring closures
      visited.set(centerIdx);
      continue;
    }

    const bool isCCW = center->getPerturbationOrder(currentPertOrder) %
                       2;  // if odd => counterclockwise @Greg Landrum?
    if (isCCW != (chiralTag == Atom::CHI_TETRAHEDRAL_CCW)) {
      // center is in wrong order => we need to inverse the offsets
      for (const auto dependentIdx : dependentIdxs) {
        zmat->invertTorsionDependence(dependentIdx);
      }
    }
    visited.set(centerIdx);
  }
}

// =============

}  // namespace DGeomHelpers
}  // namespace RDKit
