//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_REDUCEDGRAPHS_H_
#define _RD_REDUCEDGRAPHS_H_

#include <vector>
#include <cstdint>
#include <boost/dynamic_bitset.hpp>
#include <Numerics/Vector.h>

namespace RDKit {
class ROMol;

namespace ReducedGraphs {
//! \brief Generates a reduced graph representation of a molecule
/*!

  \param mol:          the molecule to be fingerprinted

  \return a new molecule

  <b>Notes:</b>
  - the caller is responsible for <tt>delete</tt>ing the result

*/
RDKIT_REDUCEDGRAPHS_EXPORT ROMol *generateMolExtendedReducedGraph(
    const ROMol &mol,
    std::vector<boost::dynamic_bitset<>> *atomTypes = nullptr);
//! \brief Generates a ErG fingerprint vector for a molecule that's already a
// reduced graph
/*!

  \param mol:           the molecule to be fingerprinted
  \param atomTypes:     [optional] contains bit vectors indicating whether each
  atom in
                        the molecule matches each type.
  \param fuzzIncrement: amount to be added to neighboring bins
  \param minPath:       minimum distance (in bonds) to be considered
  \param maxPath:       maximum distance (in bonds) to be considered

  \return the fingerprint, as a DoubleVector

  <b>Notes:</b>
  - the caller is responsible for <tt>delete</tt>ing the result

*/
RDKIT_REDUCEDGRAPHS_EXPORT RDNumeric::DoubleVector *
generateErGFingerprintForReducedGraph(
    const ROMol &mol, std::vector<boost::dynamic_bitset<>> *atomTypes = nullptr,
    double fuzzIncrement = 0.3, unsigned int minPath = 1,
    unsigned int maxPath = 15);

//! \brief Generates a ErG fingerprint vector for a molecule
/*!

  \param mol:           the molecule to be fingerprinted
  \param atomTypes:     [optional] contains bit vectors indicating whether each
  atom in
                        the molecule matches each type.
  \param fuzzIncrement: amount to be added to neighboring bins
  \param minPath:       minimum distance (in bonds) to be considered
  \param maxPath:       maximum distance (in bonds) to be considered

  \return the fingerprint, as a DoubleVector

  <b>Notes:</b>
  - the caller is responsible for <tt>delete</tt>ing the result

*/
RDKIT_REDUCEDGRAPHS_EXPORT RDNumeric::DoubleVector *getErGFingerprint(
    const ROMol &mol, std::vector<boost::dynamic_bitset<>> *atomTypes = nullptr,
    double fuzzIncrement = 0.3, unsigned int minPath = 1,
    unsigned int maxPath = 15);
}  // namespace ReducedGraphs
}  // namespace RDKit

#endif
