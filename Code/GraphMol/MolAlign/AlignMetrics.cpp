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
#include <math.h>
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

namespace RDKit {
namespace MolAlign {

double getConformerRMS(ROMol &mol, unsigned int confId1, unsigned int confId2,
                       const std::vector<unsigned int> *atomIds,
                       bool prealigned) {
  if (!prealigned) {
    static const unsigned int confIds[] = {confId1, confId2};
    const std::vector<unsigned int> confIds2(
        confIds, confIds + sizeof(confIds) / sizeof(confIds[0]));
    MolAlign::alignMolConformers(mol, atomIds, &confIds2);
  }

  const Conformer conf1 = mol.getConformer(confId1);
  const Conformer conf2 = mol.getConformer(confId2);

  const size_t nAtoms = conf1.getNumAtoms();
  double ssr = 0.0;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    RDGeom::Point3D v = conf1.getAtomPos(i) - conf2.getAtomPos(i);
    ssr += v.dotProduct(v);
  }
  ssr /= nAtoms;
  return sqrt(ssr);
}

double getBestRMS(const ROMol &refMol, ROMol &prbMol, int refCid, int prbCid,
                  std::vector<MatchVectType> *atomMaps) {
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;

  if (!atomMaps || atomMaps->size() == 0) {
    // we have to figure out all mappings between the two molecule
    const bool uniquify = true;
    const bool recursionPossible = true;
    const bool useChirality = false;
    const bool useQueryQueryMatches = true;
    unsigned int nmatches =
        SubstructMatch(refMol, prbMol, *atomMaps, uniquify, recursionPossible,
                       useChirality, useQueryQueryMatches);
    if (nmatches == 0) {
      throw MolAlignException(
          "No sub-structure match found between the probe and query mol");
    }
  }

  double bestRMS = 10000.0;
  MatchVectType bestMap;
  bool bestAtLast = false;
  BOOST_FOREACH (const MatchVectType &atomMap, *atomMaps) {
    double rms = alignMol(prbMol, refMol, prbCid, refCid, &atomMap);
    bestAtLast = false;
    if (rms < bestRMS) {
      bestRMS = rms;
      bestMap = atomMap;
      bestAtLast = true;
    }
  }
  if (!bestAtLast) {
    alignMol(prbMol, refMol, prbCid, refCid, &bestMap);
  }
  return bestRMS;
}
}
}

/*


def GetConformerRMSMatrix(mol, atomIds=None, prealigned=False):
  """ Returns the RMS matrix of the conformers of a molecule.
  As a side-effect, the conformers will be aligned to the first
  conformer (i.e. the reference) and will left in the aligned state.

  Arguments:
    - mol:     the molecule
    - atomIds: (optional) list of atom ids to use a points for
               alingment - defaults to all atoms
    - prealigned: (optional) by default the conformers are assumed
                  be unaligned and will therefore be aligned to the
                  first conformer

  Note that the returned RMS matrix is symmetrically, i.e. it is the
  lower half of the matrix, e.g. for 5 conformers:
  rmsmatrix = [ a,
                b, c,
                d, e, f,
                g, h, i, j]
  This way it can be directly used as distance matrix in e.g. Butina
  clustering.

  """
  # if necessary, align the conformers
  # Note: the reference conformer is always the first one
  rmsvals = []
  if not prealigned:
    if atomIds:
      AlignMolConformers(mol, atomIds=atomIds, RMSlist=rmsvals)
    else:
      AlignMolConformers(mol, RMSlist=rmsvals)
  else:  # already prealigned
    for i in range(1, mol.GetNumConformers()):
      rmsvals.append(GetConformerRMS(mol, 0, i, atomIds=atomIds,
prealigned=prealigned))
  # loop over the conformations (except the reference one)
  cmat = []
  for i in range(1, mol.GetNumConformers()):
    cmat.append(rmsvals[i - 1])
    for j in range(1, i):
      cmat.append(GetConformerRMS(mol, i, j, atomIds=atomIds, prealigned=True))
  return cmat
*/
