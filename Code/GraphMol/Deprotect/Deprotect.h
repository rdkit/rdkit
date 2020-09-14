//
//  Copyright (C) 2020 Brian P Kelley, Joann Prescott-Roy
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDK_DEPROTECT_LIBRARY
#define RDK_DEPROTECT_LIBRARY

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <string>
#include <memory>

namespace RDKit {
namespace Deprotect {
/*! Data for Deprotecting molecules

 Deprotects are described as reactions that remove the protecting
  group and leave behind the group being protected.

  Each DeprotectData has the following attributes:

    - <b>deprotection_class</b> functional group being protected (i.e. amine,
 alcohol, ...)
    - <b>reaction_smarts</b> the reaction smarts pattern for removing the
 protecting group
    - <b>abbreviation</b> common abbreviation for the protecting group (Boc,
 Fmoc)
    - <b>full_name</b> full name for the protecting group
    - <b> rxn </b> the reaction itself.
*/

struct RDKIT_DEPROTECT_EXPORT DeprotectData {
  std::string deprotection_class;
  std::string reaction_smarts;
  std::string abbreviation;
  std::string full_name;
  std::string example;

  std::shared_ptr<ChemicalReaction>
      rxn;  // so much easier than unique_ptr, sigh...

  DeprotectData(const std::string &deprotection_class,
                const std::string &reaction_smarts,
                const std::string &abbrevition,
		const std::string &full_name,
		const std::string &example="");

  bool operator==(const DeprotectData &other) const {
    return (deprotection_class == other.deprotection_class &&
            full_name == other.full_name &&
            abbreviation == other.abbreviation &&
            reaction_smarts == other.reaction_smarts &&
            isValid() == other.isValid());
  }

  //! Returns true if the deprotection is valid
  bool isValid() const {
    return rxn.get() != nullptr && rxn->getNumProductTemplates() == 1;
  }
};

//! Retrieves the built in list of common deprotections
RDKIT_DEPROTECT_EXPORT const std::vector<DeprotectData> &getDeprotections();

//! Deprotect a molecule
/*!
     The resulting molecule is annotated with the deprotections used (property
   DEPROTECTIONS) and the number of deprotections applied (property
   DEPROTECTIION_COUNT)

     \param mol the molecule to deprotect
     \param deprotections - a vector of deprotections to use, defaults to the
   built in deprotections.

     \return The deprotected form of the input molecule
*/
RDKIT_DEPROTECT_EXPORT std::unique_ptr<ROMol> deprotect(
    const ROMol &mol,
    const std::vector<DeprotectData> &deprotections = getDeprotections());
}  // namspace Deprotect 
} // namespace RDKit
#endif
