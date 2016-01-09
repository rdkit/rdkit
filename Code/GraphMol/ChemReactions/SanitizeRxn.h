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
#ifndef RDKIT_SANITIZERXN_H
#define RDKIT_SANITIZERXN_H

#include "Reaction.h"
#include <string>
#include <exception>

namespace RDKit {

//! class for flagging sanitization errors
class RxnSanitizeException : public std::exception {
 public:
  RxnSanitizeException(const char *msg) : _msg(msg){};
  RxnSanitizeException(const std::string &msg) : _msg(msg){};
  const char *message() const { return _msg.c_str(); };
  ~RxnSanitizeException() throw(){};

 private:
  std::string _msg;
};


namespace RxnOps {
//! If atom maps are not defined, attempt to deduce them from the RGroup
//!  labels.
void fixAtomMaps(ChemicalReaction &rxn);

//! Any dummy atom with a map but no RGroup label, should be an RGroup
//!  in RDKit's view of a reaction.
void fixRGroups(ChemicalReaction &rxn);

//! Set aromaticity on templates if possible
void fixTemplateAromaticity(ChemicalReaction &rxn);

//! merge query Hs if appropriate
void fixHs(ChemicalReaction &rxn);

typedef enum {
  SANITIZE_NONE = 0x0,
  SANITIZE_ATOM_MAPS = 0x1,
  SANITIZE_RGROUP_NAMES = 0x2,
  SANITIZE_REAGENT_AROMATICITY = 0x4,
  SANITIZE_MERGEHS = 0x8,
  SANITIZE_ALL = 0xFFFFFFF
} SanitizeRxnFlags;

//! \brief carries out a collection of tasks for cleaning up a reaction and
// ensuring
//! that it makes "chemical sense" in the context of RDKit reacitons
/*!
   This functions calls the following in sequence
     -# RxnOps::fixupAtomMaps()
     -# RxnOps::fixRGroups()
     -# MolOps::fixupTemplateAromaticity()

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
       HOWEVER, if any flag is returned in operationsPerformed, take the
       the reaction may still be suspect to its validity.
    - Aromaticity can be tricky when starting with Kekule structures that
      have query features, aromaticity works well for non-query rings, however
      certain structures (substitutions on Kekule rings that should really be
   aromatic)
      may not have enough information.
*/

void sanitizeRxn(ChemicalReaction &rxn,
                 unsigned int &operationsThatFailed,
                 unsigned int sanitizeOps = SANITIZE_ALL);
//! \overload
void sanitizeRxn(ChemicalReaction &rxn);
}
}

#endif
