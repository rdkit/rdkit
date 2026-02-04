//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <mutex>
#include <cmath>

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

std::mutex mtx;

namespace RDKit {
namespace GaussianShape {

// Bondi radii
//  can find more of these in Table 12 of this publication:
//   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
// The dummy atom radius (atomic number 0) is set to
// 2.16 in ShapeInputOptions and may be varied there, as
// may all the other radii if required, including the
// addition of atoms not covered here.
const std::map<unsigned int, double> vdw_radii = {
    {1, 1.10},   // H
    {2, 1.40},   // He
    {3, 1.81},   // Li
    {4, 1.53},   // Be
    {5, 1.92},   // B
    {6, 1.70},   // C
    {7, 1.55},   // N
    {8, 1.52},   // O
    {9, 1.47},   // F
    {10, 1.54},  // Ne
    {11, 2.27},  // Na
    {12, 1.73},  // Mg
    {13, 1.84},  // Al
    {14, 2.10},  // Si
    {15, 1.80},  // P
    {16, 1.80},  // S
    {17, 1.75},  // Cl
    {18, 1.88},  // Ar
    {19, 2.75},  // K
    {20, 2.31},  // Ca
    {31, 1.87},  // Ga
    {32, 2.11},  // Ge
    {33, 1.85},  // As
    {34, 1.90},  // Se
    {35, 1.83},  // Br
    {36, 2.02},  // Kr
    {37, 3.03},  // Rb
    {38, 2.49},  // Sr
    {49, 1.93},  // In
    {50, 2.17},  // Sn
    {51, 2.06},  // Sb
    {52, 2.06},  // Te
    {53, 1.98},  // I
    {54, 2.16},  // Xe
    {55, 3.43},  // Cs
    {56, 2.68},  // Ba
    {81, 1.96},  // Tl
    {82, 2.02},  // Pb
    {83, 2.07},  // Bi
    {84, 1.97},  // Po
    {85, 2.02},  // At
    {86, 2.20},  // Rn
    {87, 3.48},  // Fr
    {88, 2.83},  // Ra
};
constexpr double radius_color =
    1.08265;  // same radius for all feature/color "atoms", as used by the
              // PubChem code.

ShapeInput::ShapeInput(const ROMol &mol, int confId,
                       const ShapeOverlayOptions &overlayOpts) {
  extractAtoms(mol, confId);
  if (overlayOpts.optimMode != OptimMode::SHAPE_ONLY) {
    extractFeatures(mol, confId, overlayOpts);
  }
  calcNormalization(mol, confId);
  calcExtremes();
  d_selfOverlapVol =
      calcVolAndGrads(d_coords.data(), d_numAtoms, d_coords.data(), d_numAtoms,
                      overlayOpts.all_carbon_radii);
  d_selfOverlapColor = calcVolAndGrads(
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      nullptr, nullptr);
}

ShapeInput::ShapeInput(const ShapeInput &other)
    : d_coords(other.d_coords),
      d_types(other.d_types),
      d_numAtoms(other.d_numAtoms),
      d_numFeats(other.d_numFeats),
      d_selfOverlapVol(other.d_selfOverlapVol),
      d_selfOverlapColor(other.d_selfOverlapColor),
      d_extreme_points(other.d_extreme_points),
      d_normalized(other.d_normalized),
      d_canonRot(new std::array<double, 9>(*other.d_canonRot)),
      d_centroid(new std::array<double, 3>(*other.d_centroid)) {}

ShapeInput &ShapeInput::operator=(const ShapeInput &other) {
  if (this == &other) {
    return *this;
  }
  d_coords = other.d_coords;
  d_types = other.d_types;
  d_numAtoms = other.d_numAtoms;
  d_numFeats = other.d_numFeats;
  d_selfOverlapVol = other.d_selfOverlapVol;
  d_selfOverlapColor = other.d_selfOverlapColor;
  d_extreme_points = other.d_extreme_points;
  d_normalized = other.d_normalized;
  d_canonRot.reset(new std::array<double, 9>(*other.d_canonRot));
  d_centroid.reset(new std::array<double, 3>(*other.d_centroid));
  return *this;
}

void ShapeInput::normalizeCoords() {
  // For consistency, create a dummy molecule and use the same code
  // as for shapes created from a molecule.
  RWMol mol;
  Conformer *conf = new Conformer(d_coords.size() / 4);
  for (size_t i = 0, j = 0; i < d_coords.size(); i += 4, ++j) {
    Atom *atom = new Atom(6);
    mol.addAtom(atom, true, true);
    RDGeom::Point3D pos{d_coords[i], d_coords[i + 1], d_coords[i + 2]};
    conf->setAtomPos(j, pos);
  }
  auto confId = mol.addConformer(conf);
  calcNormalization(mol, confId);
  RDGeom::Transform3D canonRot;
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      canonRot.setValUnchecked(i, j, (*d_canonRot)[k]);
    }
  }
  RDGeom::Point3D trans{(*d_centroid)[0], (*d_centroid)[1], (*d_centroid)[2]};
  canonRot.TransformPoint(trans);
  canonRot.SetTranslation(trans);

  applyTransformToShape(d_coords, canonRot);
  d_normalized = true;
  // Recalculate the extremes now we've changed the coordinates.
  calcExtremes();
}

void ShapeInput::transformCoords(RDGeom::Transform3D &xform) {
  applyTransformToShape(d_coords, xform);
}

void ShapeInput::extractAtoms(const ROMol &mol, int confId) {
  d_coords.reserve(mol.getNumAtoms() * 4);
  auto conf = mol.getConformer(confId);
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() > 1) {
      unsigned int idx = atom->getIdx();
      auto &pos = conf.getAtomPos(idx);
      d_coords.push_back(pos.x);
      d_coords.push_back(pos.y);
      d_coords.push_back(pos.z);
      if (auto rad =
              vdw_radii.find(static_cast<unsigned int>(atom->getAtomicNum()));
          rad != vdw_radii.end()) {
        d_coords.push_back(KAPPA / (rad->second * rad->second));
      } else {
        throw ValueErrorException("No VdW radius for atom with Z=" +
                                  std::to_string(atom->getAtomicNum()));
      }
    }
  }
  d_numAtoms = d_coords.size() / 4;
  d_types.resize(d_numAtoms);
  d_numFeats = 0;
}

// Extract the features for the color scores, using RDKit pphore features
// for now.  Other options to be added later.
void ShapeInput::extractFeatures(const ROMol &mol, int confId,
                                 const ShapeOverlayOptions &shapeOpts) {
  const auto &pattVects = shapeOpts.d_ph4Patterns;
  unsigned pattIdx = 1;
  for (const auto &patts : pattVects) {
    for (const auto &patt : patts) {
      std::vector<MatchVectType> matches;
      {
        // recursive queries aren't thread safe.
        std::unique_lock<std::mutex> lock(mtx);
        matches = SubstructMatch(mol, *patt);
      }
      for (auto match : matches) {
        std::vector<unsigned int> ats;
        for (const auto &pr : match) {
          ats.push_back(pr.second);
        }
        auto featPos = computeFeaturePos(mol, confId, ats);
        d_types.push_back(pattIdx);
        d_coords.push_back(featPos.x);
        d_coords.push_back(featPos.y);
        d_coords.push_back(featPos.z);
        d_coords.push_back(KAPPA / (radius_color * radius_color));
        d_numFeats++;
      }
    }
    ++pattIdx;
  }
}

void ShapeInput::calcNormalization(const ROMol &mol, int confId) {
  std::unique_ptr<RDGeom::Transform3D> canonXform(
      MolTransforms::computeCanonicalTransform(mol.getConformer(confId)));
  d_canonRot.reset(
      new std::array<double, 9>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      (*d_canonRot)[k] = canonXform->getValUnchecked(i, j);
    }
  }
  d_centroid.reset(new std::array<double, 3>{0.0, 0.0, 0.0});
  for (int i = 0; i < 4 * d_numAtoms; i += 4) {
    (*d_centroid)[0] -= d_coords[i];
    (*d_centroid)[1] -= d_coords[i + 1];
    (*d_centroid)[2] -= d_coords[i + 2];
  }
  (*d_centroid)[0] /= d_numAtoms;
  (*d_centroid)[1] /= d_numAtoms;
  (*d_centroid)[2] /= d_numAtoms;
}

void ShapeInput::calcExtremes() {
  d_extreme_points = std::array<size_t, 6>{0, 0, 0, 0, 0, 0};
  for (size_t i = 0, j = 0; i < d_coords.size(); i += 4, ++j) {
    if (d_coords[i] < d_coords[4 * d_extreme_points[0]]) {
      d_extreme_points[0] = j;
    }
    if (d_coords[i] > d_coords[4 * d_extreme_points[3]]) {
      d_extreme_points[3] = j;
    }

    if (d_coords[i + 1] < d_coords[4 * d_extreme_points[1] + 1]) {
      d_extreme_points[1] = j;
    }
    if (d_coords[i + 1] > d_coords[4 * d_extreme_points[4] + 1]) {
      d_extreme_points[4] = j;
    }

    if (d_coords[i + 2] < d_coords[4 * d_extreme_points[2] + 2]) {
      d_extreme_points[2] = j;
    }
    if (d_coords[i + 2] > d_coords[4 * d_extreme_points[5] + 2]) {
      d_extreme_points[5] = j;
    }
  }
}

RDGeom::Point3D computeFeaturePos(const ROMol &mol, int confId,
                                  const std::vector<unsigned int> &ats) {
  RDGeom::Point3D featPos;
  auto &conf = mol.getConformer(confId);
  for (const auto at : ats) {
    featPos += conf.getAtomPos(at);
  }
  featPos /= ats.size();
  return featPos;
}

void writeCoords(const std::vector<DTYPE> &shape, const std::string &label,
                 char lineEnd) {
  std::cout << label << " :: ";
  if (lineEnd == '\n') {
    std::cout << lineEnd;
  }
  for (size_t i = 0; i < shape.size(); i += 4) {
    std::cout << shape[i] << "," << shape[i + 1] << "," << shape[i + 2]
              << lineEnd;
  }
}

void writeCoords(const DTYPE *shape, unsigned int numPts,
                 const std::string &label, char lineEnd) {
  std::cout << label << " :: ";
  if (lineEnd == '\n') {
    std::cout << lineEnd;
  }
  for (unsigned int i = 0; i < numPts * 4; i += 4) {
    std::cout << shape[i] << "," << shape[i + 1] << "," << shape[i + 2]
              << lineEnd;
  }
}

void copyTransform(const RDGeom::Transform3D &src, RDGeom::Transform3D &dest) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      dest.setValUnchecked(i, j, src.getValUnchecked(i, j));
    }
  }
}

void applyTransformToShape(std::vector<DTYPE> &shape,
                           RDGeom::Transform3D &xform) {
  for (size_t i = 0; i < shape.size(); i += 4) {
    RDGeom::Point3D pos{shape[i], shape[i + 1], shape[i + 2]};
    xform.TransformPoint(pos);
    shape[i] = pos.x;
    shape[i + 1] = pos.y;
    shape[i + 2] = pos.z;
  }
}

void applyTransformToShape(const DTYPE *inShape, DTYPE *outShape,
                           size_t numPoints, RDGeom::Transform3D &xform) {
  for (size_t i = 0; i < 4 * numPoints; i += 4) {
    RDGeom::Point3D pos{inShape[i], inShape[i + 1], inShape[i + 2]};
    xform.TransformPoint(pos);
    outShape[i] = pos.x;
    outShape[i + 1] = pos.y;
    outShape[i + 2] = pos.z;
    outShape[i + 3] = inShape[i + 3];
  }
}

void translateShape(std::vector<DTYPE> &shape,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < shape.size(); i += 4) {
    shape[i] += translation.x;
    shape[i + 1] += translation.y;
    shape[i + 2] += translation.z;
  }
}

void translateShape(const DTYPE *inShape, DTYPE *outShape, size_t numPoints,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < 4 * numPoints; i += 4) {
    outShape[i] = inShape[i] + translation.x;
    outShape[i + 1] = inShape[i + 1] + translation.y;
    outShape[i + 2] = inShape[i + 2] + translation.z;
    outShape[i + 3] = inShape[i + 3];
  }
}
}  // namespace GaussianShape
}  // namespace RDKit
