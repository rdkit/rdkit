//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_ALIGNMETRICS_H_
#define _RD_ALIGNMETRICS_H_

#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <vector>

namespace RDKit {
namespace MolAlign {

//! Alignment metrics

//! Compute the RMSD between two conformations
/*!
  Returns the RMS between two conformations.
  By default, the conformers will be aligned to the first conformer
  of the molecule (i.e. the reference) before RMS calculation and,
  as a side-effect, will be left in the aligned state.

  \param mol        molecule with several conformations
  \param confId1    the id of the first conformer
  \param confId2    the id of the second conformer
  \param atomIds    list of atom ids to use as points for alignment (defaults
                    to all atoms)
  \param prealigned by default the conformers are assumed to be unaligned and
                    will therefore be aligned to the first conformer

  <b>Returns</b>
  RMSD value
*/
double getConformerRMS(ROMol &mol, unsigned int confId1, unsigned int confId2,
                       const std::vector<unsigned int> *atomIds = 0,
                       bool prealigned = false);

//! Compute the optimal RMSE between two molecules
/*!
  Returns the optimal RMS for aligning two molecules, taking
  symmetry into account. As a side-effect, the probe molecule is
  left in the aligned state.

  \param refMol       the reference molecule
  \param prbMol       the molecule to be aligned to the reference
  \param refCid       (optional) reference conformation to use
  \param prbCid       (optional) probe conformation to use
  \param atomMaps     (optional) a list of atomMaps between the two molecules.
                      If not provided, these will be generated using a
                      substructure search.

  Note:
  This function will attempt to align all permutations of matching atom
  orders in both molecules, for some molecules it will lead to 'combinatorial
  explosion' especially if hydrogens are present.
  Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing the
  atom order.
  """
*/
double getBestRMS(const ROMol &refMol, ROMol &prbMol, int probeConfId = -1,
                  int refConfId = -1, std::vector<MatchVectType> *atomMaps = 0);
}
}
#endif
