//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_SANITIZERXN_H
#define RDKIT_SANITIZERXN_H

#include "Reaction.h"
#include <GraphMol/MolOps.h>
#include <string>
#include <exception>

namespace RDKit {

//! class for flagging sanitization errors
class RDKIT_CHEMREACTIONS_EXPORT RxnSanitizeException : public std::exception {
 public:
  RxnSanitizeException(const char *msg) : _msg(msg){};
  RxnSanitizeException(const std::string &msg) : _msg(msg){};
  const char *what() const noexcept override { return _msg.c_str(); };
  ~RxnSanitizeException() noexcept {};

 private:
  std::string _msg;
};

namespace RxnOps {
//! Any dummy atom with a map but no RGroup label, should be an RGroup
//!  in RDKit's view of a reaction.
//!  See if these atoms can be salvaged into RGroups.
RDKIT_CHEMREACTIONS_EXPORT void fixRGroups(ChemicalReaction &rxn);

//! If atom maps are not defined on rgroups, attempt to deduce them from the
//! RGroup
//!  labels, or add new ones if possible.
RDKIT_CHEMREACTIONS_EXPORT void fixAtomMaps(ChemicalReaction &rxn);

//! Adjusts the reactant templates to properly match reagents
RDKIT_CHEMREACTIONS_EXPORT void adjustTemplates(
    ChemicalReaction &rxn, const MolOps::AdjustQueryParameters &params);

//! merge query Hs if appropriate
RDKIT_CHEMREACTIONS_EXPORT void fixHs(ChemicalReaction &rxn);

// Default adjustment parameters for matching reagents
inline const MolOps::AdjustQueryParameters DefaultRxnAdjustParams() {
  MolOps::AdjustQueryParameters params;
  params.adjustDegree = false;
  params.adjustDegreeFlags = MolOps::ADJUST_IGNOREALL;
  params.adjustRingCount = false;
  params.adjustRingCountFlags = MolOps::ADJUST_IGNOREALL;
  params.makeDummiesQueries = false;
  params.aromatizeIfPossible = true;
  return params;
}

// Default adjustment parameters for ChemDraw style matching of reagents
//  -- deprecated - renamed MatchOnlyAtRgroupsAdjustParams
//  -- this doesn't match sciquest style searching
inline const MolOps::AdjustQueryParameters ChemDrawRxnAdjustParams() {
  BOOST_LOG(rdWarningLog)
      << " deprecated -- please use MatchOnlyAtRgroupsAdjustParams instead"
      << std::endl;
  MolOps::AdjustQueryParameters params;
  params.adjustDegree = true;
  params.adjustDegreeFlags = MolOps::ADJUST_IGNOREDUMMIES;
  params.adjustRingCount = false;
  params.adjustRingCountFlags = MolOps::ADJUST_IGNORENONE;
  params.makeDummiesQueries = false;
  params.aromatizeIfPossible = true;
  return params;
}

inline const MolOps::AdjustQueryParameters MatchOnlyAtRgroupsAdjustParams() {
  MolOps::AdjustQueryParameters params;
  params.adjustDegree = true;
  params.adjustDegreeFlags = MolOps::ADJUST_IGNOREDUMMIES;
  params.adjustRingCount = false;
  params.adjustRingCountFlags = MolOps::ADJUST_IGNORENONE;
  params.makeDummiesQueries = false;
  params.aromatizeIfPossible = true;
  return params;
}

typedef enum {
  SANITIZE_NONE = 0x0,
  SANITIZE_RGROUP_NAMES = 0x1,
  SANITIZE_ATOM_MAPS = 0x2,
  SANITIZE_ADJUST_REACTANTS = 0x4,
  SANITIZE_MERGEHS = 0x8,
  SANITIZE_ALL = 0xFFFFFFFF
} SanitizeRxnFlags;

//! \brief carries out a collection of tasks for cleaning up a reaction and
// ensuring
//! that it makes "chemical sense" in the context of RDKit reacitons
/*!
   This functions calls the following in sequence
     -# RxnOps::fixRGroups()
     -# RxnOps::fixupAtomMaps()
     -# RxnOps::fixupTemplateAromaticity()
     -# RxnOps::mergeHs()

   \param rxn : the ChemicalReaction to be cleaned

   \param operationThatFailed : the first (if any) sanitization operation that
   fails is set here.
                                The values are taken from the \c SanitizeFlags
   enum.
                                On success, the value is  \c
   SanitizeFlags::SANITIZE_NONE

   \param sanitizeOps : the bits here are used to set which sanitization
   operations are carried
                        out. The elements of the \c SanitizeFlags enum define
   the operations.

   <b>Notes:</b>
    - This attempts to fix known issues with certain reaction drawers.
       HOWEVER, if any flag is returned in operationsPerformed,
       the reaction may still be suspect to its validity.
    - Aromaticity can be tricky when starting with Kekule structures that
      have query features, aromaticity works well for non-query rings, however
      certain structures (substitutions on Kekule rings that should really be
      aromatic) may not have enough information.
*/

RDKIT_CHEMREACTIONS_EXPORT void sanitizeRxn(
    ChemicalReaction &rxn, unsigned int &operationsThatFailed,
    unsigned int sanitizeOps = SANITIZE_ALL,
    const MolOps::AdjustQueryParameters &params = DefaultRxnAdjustParams());
//! \overload
RDKIT_CHEMREACTIONS_EXPORT void sanitizeRxn(
    ChemicalReaction &rxn,
    const MolOps::AdjustQueryParameters &params = DefaultRxnAdjustParams());

}  // namespace RxnOps
}  // namespace RDKit

#endif
