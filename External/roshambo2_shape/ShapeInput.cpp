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

#include <cmath>
#include <mutex>

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "ShapeInput.h"
#include "roshambo2/roshambo2/backends/cpp_src/cpp_helper_functions.h"

#include <boost/locale/boundary/types.hpp>

std::mutex mtx;

namespace RDKit {
namespace ShapeAlign {

ShapeInput::ShapeInput(const ROMol &mol, int confId,
                       const ShapeOverlayOptions &overlayOpts) {
  extractAtoms(mol, confId);
  if (overlayOpts.d_useColors) {
    extractFeatures(mol, confId, overlayOpts);
  }
  calcNormalization(mol, confId);
  d_selfOverlapVol =
      volume(d_coords.data(), d_numAtoms, d_coords.data(), d_numAtoms);
  d_selfOverlapColor = volume_color(
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      overlayOpts.d_rMat.data(), overlayOpts.d_pMat.data(),
      overlayOpts.d_nTypes);
}

ShapeInput::ShapeInput(const ShapeInput &other)
    : d_coords(other.d_coords),
      d_types(other.d_types),
      d_numAtoms(other.d_numAtoms),
      d_numFeats(other.d_numFeats),
      d_selfOverlapVol(other.d_selfOverlapVol),
      d_selfOverlapColor(other.d_selfOverlapColor),
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
  d_normalized = other.d_normalized;
  d_canonRot.reset(new std::array<double, 9>(*other.d_canonRot));
  d_centroid.reset(new std::array<double, 3>(*other.d_centroid));
  return *this;
}

void ShapeInput::normalizeCoords() {
  if (!d_canonRot) {
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
  }
  RDGeom::Transform3D canonRot;
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      canonRot.setValUnchecked(i, j, (*d_canonRot)[k]);
    }
  }
  RDGeom::Point3D trans{(*d_centroid)[0], (*d_centroid)[1], (*d_centroid)[2]};
  canonRot.TransformPoint(trans);
  canonRot.SetTranslation(trans);

  for (size_t i = 0; i < d_coords.size(); i += 4) {
    RDGeom::Point3D pos{d_coords[i], d_coords[i + 1], d_coords[i + 2]};
    canonRot.TransformPoint(pos);
    d_coords[i] = pos.x;
    d_coords[i + 1] = pos.y;
    d_coords[i + 2] = pos.z;
  }
  d_normalized = true;
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
      d_coords.push_back(1.0);
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
        d_coords.push_back(1.0);
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

RDGeom::Transform3D quatTransToTransform(const DTYPE *quat,
                                         const DTYPE *trans) {
  RDGeom::Transform3D xform;
  xform.setToIdentity();
  DTYPE rot[3][3];
  std::array<DTYPE, 4> tquat{quat[0], quat[1], quat[2], quat[3]};
  // Use Roshambo2's function.  Transform3D::SetRotationFromQuaternion
  // gives a different answer.
  quaternion_to_rotation_matrix(tquat, rot);
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      xform.setValUnchecked(i, j, rot[i][j]);
    }
  }
  xform.SetTranslation(RDGeom::Point3D{trans[0], trans[1], trans[2]});
  return xform;
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

}  // namespace ShapeAlign
}  // namespace RDKit
