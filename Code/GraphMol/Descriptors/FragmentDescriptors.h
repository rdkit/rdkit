//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file FragmentDescriptors.h

  \brief C++ port of the fr_* fragment count descriptors from
  rdkit.Chem.Fragments.

  The SMARTS definitions are those from $RDBASE/Data/FragmentDescriptors.csv
  and are compiled into the library. Names follow the python module, i.e.
  "=" and "-" in the names from the data file are replaced by "_".

*/
#include <RDGeneral/export.h>
#ifndef RD_FRAGMENTDESCRIPTORS_H
#define RD_FRAGMENTDESCRIPTORS_H

#include <string>
#include <vector>

namespace RDKit {
class ROMol;
namespace Descriptors {

const std::string FragmentDescriptorsVersion = "1.0.0";

//! returns the names of the fragment descriptors (fr_C_O, fr_C_O_noCOO, ...)
//! in the order used by calcFragmentDescriptors()
RDKIT_DESCRIPTORS_EXPORT const std::vector<std::string> &
getFragmentDescriptorNames();

//! returns the SMARTS definitions of the fragment descriptors, in the same
//! order as getFragmentDescriptorNames()
RDKIT_DESCRIPTORS_EXPORT const std::vector<std::string> &
getFragmentDescriptorSmarts();

//! returns the one-line descriptions of the fragment descriptors, in the
//! same order as getFragmentDescriptorNames()
RDKIT_DESCRIPTORS_EXPORT const std::vector<std::string> &
getFragmentDescriptorDescriptions();

//! returns the number of matches of each fragment pattern
/*!
  \param mol          the molecule of interest
  \param countUnique  if set (the default), only unique matches are counted
*/
RDKIT_DESCRIPTORS_EXPORT std::vector<unsigned int> calcFragmentDescriptors(
    const ROMol &mol, bool countUnique = true);

//! returns the number of matches of a single fragment pattern
/*!
  \param mol          the molecule of interest
  \param name         the descriptor name, e.g. "fr_benzene"
  \param countUnique  if set (the default), only unique matches are counted

  throws a KeyErrorException if \c name is not a known fragment descriptor
*/
RDKIT_DESCRIPTORS_EXPORT unsigned int calcFragmentDescriptor(
    const ROMol &mol, const std::string &name, bool countUnique = true);

}  // namespace Descriptors
}  // namespace RDKit
#endif
