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
#ifndef ENUMERATETYPES_H
#define ENUMERATETYPES_H

#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace EnumerationTypes {
//! BBS - Helper typedef for holding building blocks for reactions
//!   holds vectors of reagents for each reactant in a Reaction
typedef std::vector<MOL_SPTR_VECT> BBS;

//! RGROUPS Helper typedef for indexing into the BBS vector
//!  - The indices into the BBS molecule list to create a product
//!  Example
//!   RGROUPS groups;
//!   groups.push_back(10);
//!   groups.push_back(5);
//!
//!   Will create a product from the following building blocks:
//!   MOL_SPTR_VECT building_blocks;
//!   building_blocks.push_back( BBS[0][groups[0] );
//!   building_blocks.push_back( BBS[1][groups[1] );
//!    rxn.runReactants( building_blocks );
typedef std::vector<boost::uint64_t> RGROUPS;
}  // namespace EnumerationTypes
}  // namespace RDKit
#endif
