//
//  Copyright (C) 2017-2026 Sereina Riniker and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_TORSIONPREFERENCES_H_
#define _RD_TORSIONPREFERENCES_H_
#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <tuple>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>

namespace RDKit {
class ROMol;
}  // namespace RDKit

namespace ForceFields {
namespace CrystalFF {

//! A structure used to the experimental torsion patterns
struct RDKIT_FORCEFIELDHELPERS_EXPORT ExpTorsionAngle {
  unsigned int torsionIdx;
  std::string smarts;
  std::vector<double> V;
  std::vector<int> signs;
  std::unique_ptr<const RDKit::ROMol> dp_pattern;
  unsigned int idx[4];
};

struct RDKIT_FORCEFIELDHELPERS_EXPORT GaussianExpTorsionAngle {
  std::size_t torsionIdx;
  std::string smarts;
  std::vector<double> heights;
  std::vector<double> positions;
  std::vector<double> widths;
  std::unique_ptr<const RDKit::ROMol> dp_pattern;
  unsigned int idx[4];
};

namespace ETKDGForceConsts {
struct Params {
  double distance{1.0};
  double fourthDim{1.0};
  double chiral{1.0};
  double kTermAngle{1.0};
  double kTermImproper{10.0};
  double kTermTorsion{100.0};
  double etTermScaling{1.0};
};

namespace SEQ {
constexpr Params Cosine;
constexpr Params Gaussian;
}  // namespace SEQ
namespace AIO {
constexpr Params Cosine = {.distance = 2.15,
                           .fourthDim = 2.15,
                           .kTermAngle = 0.1,
                           .kTermImproper = .001,
                           .kTermTorsion = 2.15,
                           .etTermScaling = 0.05};

constexpr Params Gaussian = {.distance = 2.15,
                             .fourthDim = 2.15,
                             .kTermAngle = 0.1,
                             .kTermImproper = .001,
                             .kTermTorsion = .25,
                             .etTermScaling = 0.05};

}  // namespace AIO
}  // namespace ETKDGForceConsts

using CosineExp_T = std::pair<std::vector<int>, std::vector<double>>;
using GaussianExp_T = std::tuple<std::vector<double>, std::vector<double>,
                                 std::vector<double>, double>;

using TorsionLookup = std::vector<std::vector<double>>;
template <typename T>
concept TorsionAngleType = std::is_same_v<T, ExpTorsionAngle> ||
    std::is_same_v<T, GaussianExpTorsionAngle>;

template <typename T>
concept TorsionParamType =
    std::is_same_v<T, CosineExp_T> || std::is_same_v<T, GaussianExp_T>;

template <TorsionParamType T>
using MappedAngle_T =
    std::conditional_t<std::is_same_v<T, CosineExp_T>, ExpTorsionAngle,
                       GaussianExpTorsionAngle>;

template <TorsionParamType T = CosineExp_T>
struct CrystalFFDetails {
  std::vector<std::vector<int>> expTorsionAtoms;
  std::vector<T> expTorsionAngles;
  std::vector<std::size_t> torsionIdx;
  std::vector<std::vector<int>> improperAtoms;
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::vector<int>> angles;
  std::vector<int> atomNums;
  double boundsMatForceScaling;
  boost::dynamic_bitset<> constrainedAtoms;
  double *distMat;
  ETKDGForceConsts::Params forceConsts;
  std::vector<RDKit::DGeomHelpers::Path14Configuration> path14Configs;
  TorsionLookup cosPhiToEnergy;
  TorsionLookup cosPhiToGrad;
};

//! Get the experimental torsional angles in a molecule
template <TorsionParamType T>
RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &mol, CrystalFFDetails<T> &details,
    bool useExpTorsions = false, bool useSmallRingTorsions = false,
    bool useMacrocycleTorsions = false, bool useBasicKnowledge = false,
    unsigned int version = 2, bool verbose = false);

//! \overload
template <TorsionParamType T>
RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &mol, CrystalFFDetails<T> &details,
    std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                           const MappedAngle_T<T> *>> &torsionBonds,
    bool useExpTorsions = false, bool useSmallRingTorsions = false,
    bool useMacrocycleTorsions = false, bool useBasicKnowledge = false,
    unsigned int version = 2, bool verbose = false);

//! Populate the lookuptable for the minimizations
RDKIT_FORCEFIELDHELPERS_EXPORT void populateRefTable(
    CrystalFFDetails<GaussianExp_T> &details);

}  // namespace CrystalFF
}  // namespace ForceFields

#endif
