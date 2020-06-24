//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_ALIGNMOLECULES_H_
#define _RD_ALIGNMOLECULES_H_

#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <vector>

namespace RDKit {
typedef std::vector<std::pair<int, int>> MatchVectType;

class Conformer;
class ROMol;
namespace MolAlign {
class RDKIT_MOLALIGN_EXPORT MolAlignException : public std::exception {
 public:
  //! construct with an error message
  MolAlignException(const char *msg) : _msg(msg){};
  //! construct with an error message
  MolAlignException(const std::string msg) : _msg(msg){};
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); };
  ~MolAlignException() noexcept {};

 private:
  std::string _msg;
};

//! Alignment functions

//! Compute the transformation required to align a molecule
/*!
  The 3D transformation required to align the specied conformation in the probe
  molecule
  to a specified conformation in the reference molecule is computed so that the
  root mean
  squared distance between a specified set of atoms is minimized

  \param prbMol    molecule that is to be aligned
  \param refMol    molecule used as the reference for the alignment
  \param trans     storage for the computed transform
  \param prbCid    ID of the conformation in the probe to be used
                   for the alignment (defaults to first conformation)
  \param refCid    ID of the conformation in the ref molecule to which
                   the alignment is computed (defaults to first conformation)
  \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                   used to compute the alignments. If this mapping is
                   not specified an attempt is made to generate on by
                   substructure matching
  \param weights   Optionally specify weights for each of the atom pairs
  \param reflect   if true reflect the conformation of the probe molecule
  \param maxIters  maximum number of iteration used in mimizing the RMSD

  <b>Returns</b>
  RMSD value
*/
RDKIT_MOLALIGN_EXPORT double getAlignmentTransform(
    const ROMol &prbMol, const ROMol &refMol, RDGeom::Transform3D &trans,
    int prbCid = -1, int refCid = -1, const MatchVectType *atomMap = nullptr,
    const RDNumeric::DoubleVector *weights = nullptr, bool reflect = false,
    unsigned int maxIters = 50);

//! Optimally (minimum RMSD) align a molecule to another molecule
/*!
  The 3D transformation required to align the specied conformation in the probe
  molecule
  to a specified conformation in the reference molecule is computed so that the
  root mean
  squared distance between a specified set of atoms is minimized. This
  transforms is them
  applied to the specified conformation in the probe molecule

  \param prbMol    molecule that is to be aligned
  \param refMol    molecule used as the reference for the alignment
  \param prbCid    ID of the conformation in the probe to be used
                   for the alignment (defaults to first conformation)
  \param refCid    ID of the conformation in the ref molecule to which
                   the alignment is computed (defaults to first conformation)
  \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                   used to compute the alignments. If this mapping is
                   not specified an attempt is made to generate on by
                   substructure matching
  \param weights   Optionally specify weights for each of the atom pairs
  \param reflect   if true reflect the conformation of the probe molecule
  \param maxIters  maximum number of iteration used in mimizing the RMSD

  <b>Returns</b>
  RMSD value
*/
RDKIT_MOLALIGN_EXPORT double alignMol(
    ROMol &prbMol, const ROMol &refMol, int prbCid = -1, int refCid = -1,
    const MatchVectType *atomMap = nullptr,
    const RDNumeric::DoubleVector *weights = nullptr, bool reflect = false,
    unsigned int maxIters = 50);

//! Returns the optimal RMS for aligning two molecules, taking
//  symmetry into account. As a side-effect, the probe molecule is
//  left in the aligned state.
/*!
  This function will attempt to align all permutations of matching atom
  orders in both molecules, for some molecules it will lead to 'combinatorial
  explosion' especially if hydrogens are present.
  Use 'RDKit::MolAlign::alignMol' to align molecules without changing the
  atom order.

  \param probeMol   the molecule to be aligned to the reference
  \param refMol     the reference molecule
  \param probeId    (optional) probe conformation to use
  \param refId      (optional) reference conformation to use
  \param map        (optional) a vector of vectors of pairs of atom IDs
                    (probe AtomId, ref AtomId) used to compute the alignments.
                    If not provided, these will be generated using a
                    substructure search.
  \param maxMatches (optional) if map is empty, this will be the max number of
                    matches found in a SubstructMatch().

  <b>Returns</b>
  Best RMSD value found
*/
RDKIT_MOLALIGN_EXPORT double getBestRMS(
    ROMol &probeMol, ROMol &refMol, int probeId = -1, int refId = -1,
    const std::vector<MatchVectType> &map = std::vector<MatchVectType>(),
    int maxMatches = 1e6);

//! Returns the RMS between two molecules, taking symmetry into account.
/*!
  This function will attempt to align all permutations of matching atom
  orders in both molecules, for some molecules it will lead to 'combinatorial
  explosion' especially if hydrogens are present.

  \param probeMol   the molecule to be aligned to the reference
  \param refMol     the reference molecule
  \param probeId    (optional) probe conformation to use
  \param refId      (optional) reference conformation to use
  \param map        (optional) a vector of vectors of pairs of atom IDs
                    (probe AtomId, ref AtomId) used to compute the alignments.
                    If not provided, these will be generated using a
                    substructure search.
  \param maxMatches (optional) if map is empty, this will be the max number of
                    matches found in a SubstructMatch().
  \param weights    (optional) weights for each pair of atom.

  <b>Returns</b>
  Best RMSD value found
*/
RDKIT_MOLALIGN_EXPORT double CalcRMS(
    ROMol &probeMol, const ROMol &refMol, int prbCid = -1, int refCid = -1,
    const std::vector<MatchVectType> &map = std::vector<MatchVectType>(),
    int maxMatches = 1e6, const RDNumeric::DoubleVector *weights = nullptr);

//! Align the conformations of a molecule using a common set of atoms. If
// the molecules contains queries, then the queries must also match exactly.

/*!
  \param mol       The molecule of interest.
  \param atomIds   vector of atoms to be used to generate the alignment.
                   All atoms will be used is not specified
  \param confIds   vector of conformations to align - defaults to all
  \param weights   vector of weights to applied to particular atom pairs
                   defaults to all weights = 1
  \param reflect   toggles reflecting (about the origin) the alignment
  \param maxIters  the maximum number of iterations to attempt
  \param RMSlist   if nonzero, this will be used to return the RMS values
                   between the reference conformation and the other aligned
                   conformations
*/
RDKIT_MOLALIGN_EXPORT void alignMolConformers(
    ROMol &mol, const std::vector<unsigned int> *atomIds = nullptr,
    const std::vector<unsigned int> *confIds = nullptr,
    const RDNumeric::DoubleVector *weights = nullptr, bool reflect = false,
    unsigned int maxIters = 50, std::vector<double> *RMSlist = nullptr);
}  // namespace MolAlign
}  // namespace RDKit
#endif
