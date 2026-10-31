//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file SAScore.h

  \brief synthetic accessibility score of Ertl and Schuffenhauer.

*/
#include <RDGeneral/export.h>
#ifndef RD_SASCORE_H
#define RD_SASCORE_H

#include <cstdint>
#include <string>
#include <vector>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! maps radius-2 Morgan bit ids onto the fragment contributions of
//! calcSAScore()
class RDKIT_DESCRIPTORS_EXPORT SAScoreFragmentTable {
 public:
  //! contribution assigned to fragments missing from the table
  static constexpr double defaultScore = -4.0;

  //! reads a table written by $RDBASE/Data/SA_Score/build_fpscores_bin.py
  /*!
    \param filename path to a binary fragment table

    Throws a BadFileException if the file cannot be opened and a
    ValueErrorException if its contents are not a valid table.
  */
  explicit SAScoreFragmentTable(const std::string &filename);

  //! contribution of a Morgan bit id, or defaultScore if it is not present
  double score(std::uint32_t bitId) const;

  //! number of fragments in the table
  std::size_t size() const { return d_keys.size(); }

 private:
  std::vector<std::uint32_t> d_keys;   //!< ascending, searched with lower_bound
  std::vector<std::uint16_t> d_index;  //!< parallel to d_keys, indexes d_scores
  std::vector<float> d_scores;
};

//! the table shipped in $RDBASE/Data/SA_Score/fpscores.bin
/*!
  Loaded on first use and shared. Requires RDBASE to be set.
*/
RDKIT_DESCRIPTORS_EXPORT const SAScoreFragmentTable &
defaultSAScoreFragmentTable();

//! calculates the synthetic accessibility score of a molecule
/*!
  Implements the method of Ertl and Schuffenhauer, J. Cheminform. 1:8 (2009),
  with the macrocycle-penalty and fingerprint-density adjustments of
  $RDBASE/Contrib/SA_Score/sascorer.py.

  Stereocenters are counted with Chirality::findPotentialStereo(), whereas
  sascorer.py uses the legacy stereo perception. The new perception finds some
  centers the legacy code misses, dependent stereocenters among them, so the
  two give different scores for molecules that have such centers.

  \param mol   the molecule of interest
  \param table the fragment contributions to use

  \return a score between 1 (easy to make) and 10 (hard to make)

  Throws a ValueErrorException if the molecule has no atoms.
*/
RDKIT_DESCRIPTORS_EXPORT double calcSAScore(const ROMol &mol,
                                            const SAScoreFragmentTable &table);
//! \overload uses defaultSAScoreFragmentTable()
RDKIT_DESCRIPTORS_EXPORT double calcSAScore(const ROMol &mol);

}  // end of namespace Descriptors
}  // end of namespace RDKit
#endif
