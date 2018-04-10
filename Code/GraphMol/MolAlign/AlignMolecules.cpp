// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignMolecules.h"
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

namespace RDKit {
namespace MolAlign {
double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid, int refCid,
                             const MatchVectType *atomMap,
                             const RDNumeric::DoubleVector *weights,
                             bool reflect, unsigned int maxIterations) {
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  const Conformer &prbCnf = prbMol.getConformer(prbCid);
  const Conformer &refCnf = refMol.getConformer(refCid);
  if (atomMap == nullptr) {
    // we have to figure out the mapping between the two molecule
    MatchVectType match;
    const bool recursionPossible = true;
    const bool useChirality = false;
    const bool useQueryQueryMatches = true;
    if (SubstructMatch(refMol, prbMol, match, recursionPossible, useChirality, useQueryQueryMatches)) {
      MatchVectType::const_iterator mi;
      for (mi = match.begin(); mi != match.end(); mi++) {
        prbPoints.push_back(&prbCnf.getAtomPos(mi->first));
        refPoints.push_back(&refCnf.getAtomPos(mi->second));
      }
    } else {
      throw MolAlignException(
          "No sub-structure match found between the probe and query mol");
    }
  } else {
    MatchVectType::const_iterator mi;
    for (mi = atomMap->begin(); mi != atomMap->end(); mi++) {
      prbPoints.push_back(&prbCnf.getAtomPos(mi->first));
      refPoints.push_back(&refCnf.getAtomPos(mi->second));
    }
  }
  double ssr = RDNumeric::Alignments::AlignPoints(
      refPoints, prbPoints, trans, weights, reflect, maxIterations);
  ssr /= (prbPoints.size());
  return sqrt(ssr);
}

double alignMol(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
                const MatchVectType *atomMap,
                const RDNumeric::DoubleVector *weights, bool reflect,
                unsigned int maxIterations) {
  RDGeom::Transform3D trans;
  double res = getAlignmentTransform(prbMol, refMol, trans, prbCid, refCid,
                                     atomMap, weights, reflect, maxIterations);
  // now transform the relevant conformation on prbMol
  Conformer &conf = prbMol.getConformer(prbCid);
  MolTransforms::transformConformer(conf, trans);
  return res;
}

double getBestRMS(ROMol &probeMol, ROMol &refMol, int probeId, int refId,
                  const std::vector<MatchVectType> &map, int maxMatches) {
  std::vector<MatchVectType> matches = map;
  if (matches.empty()) {
    bool uniquify = false;
    bool recursionPossible = true;
    bool useChirality = false;
    bool useQueryQueryMatches = false;

    SubstructMatch(refMol, probeMol, matches, uniquify, recursionPossible,
                   useChirality, useQueryQueryMatches, maxMatches);

    if (matches.empty()) {
      throw MolAlignException(
          "No sub-structure match found between the reference and probe mol");
    }

    if (matches.size() > 1e6) {
      std::string name;
      probeMol.getPropIfPresent(common_properties::_Name, name);
      std::cerr << "Warning in " << __FUNCTION__ << ": " << matches.size()
                << " matches detected for molecule " << name << ", this may "
                << "lead to a performance slowdown.\n";
    }
  }

  double bestRMS = 1.e300;
  MatchVectType &bestMatch = matches[0];
  for (auto &matche : matches) {
    double rms = alignMol(probeMol, refMol, probeId, refId, &matche);
    if (rms < bestRMS) {
      bestRMS = rms;
      bestMatch = matche;
    }
  }

  // Perform a final alignment to the best alignment...
  if (&bestMatch != &matches.back())
    alignMol(probeMol, refMol, probeId, refId, &bestMatch);
  return bestRMS;
}

void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
                        const std::vector<unsigned int> *atomIds = nullptr) {
  unsigned int na = conf.getNumAtoms();
  pts.clear();
  if (atomIds == nullptr) {
    unsigned int ai;
    pts.reserve(na);
    for (ai = 0; ai < na; ++ai) {
      pts.push_back(&conf.getAtomPos(ai));
    }
  } else {
    pts.reserve(atomIds->size());
    std::vector<unsigned int>::const_iterator cai;
    for (cai = atomIds->begin(); cai != atomIds->end(); cai++) {
      pts.push_back(&conf.getAtomPos(*cai));
    }
  }
}

void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds,
                        const std::vector<unsigned int> *confIds,
                        const RDNumeric::DoubleVector *weights, bool reflect,
                        unsigned int maxIters, std::vector<double> *RMSlist) {
  if (mol.getNumConformers() == 0) {
    // nothing to be done ;
    return;
  }

  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  int cid = -1;
  if ((confIds != nullptr) && (confIds->size() > 0)) {
    cid = confIds->front();
  }
  const Conformer &refCnf = mol.getConformer(cid);
  _fillAtomPositions(refPoints, refCnf, atomIds);

  // now loop throught the remaininf conformations and transform them
  RDGeom::Transform3D trans;
  double ssd;
  if (confIds == nullptr) {
    unsigned int i = 0;
    ROMol::ConformerIterator cnfi;
    // Conformer *conf;
    for (cnfi = mol.beginConformers(); cnfi != mol.endConformers(); cnfi++) {
      // conf = (*cnfi);
      i += 1;
      if (i == 1) {
        continue;
      }
      _fillAtomPositions(prbPoints, *(*cnfi), atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(*(*cnfi), trans);
    }
  } else {
    std::vector<unsigned int>::const_iterator cai;
    unsigned int i = 0;
    for (cai = confIds->begin(); cai != confIds->end(); cai++) {
      i += 1;
      if (i == 1) {
        continue;
      }
      Conformer &conf = mol.getConformer(*cai);
      _fillAtomPositions(prbPoints, conf, atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(conf, trans);
    }
  }
}
}
}
