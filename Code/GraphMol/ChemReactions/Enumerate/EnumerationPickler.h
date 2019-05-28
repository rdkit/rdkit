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
#ifndef ENUMERATIONPICKLER_H
#define ENUMERATIONPICKLER_H

#include "EnumerationStrategyBase.h"

namespace RDKit {
namespace EnumerationStrategyPickler {
//! pickles a EnumerationStrategy and adds the results to a stream \c ss
RDKIT_CHEMREACTIONS_EXPORT void pickle(
    const boost::shared_ptr<EnumerationStrategyBase> &enumerator,
    std::ostream &ss);
RDKIT_CHEMREACTIONS_EXPORT void pickle(
    const boost::shared_ptr<EnumerationStrategyBase> &enumerator,
    std::string &s);

//! constructs a EnumerationStrategy from a pickle stored in a string
//!  Since an EnumerationStrategyBase is polymorphic, this must return
//!  a shared pointer to the EnumerationStrategyBase
RDKIT_CHEMREACTIONS_EXPORT boost::shared_ptr<EnumerationStrategyBase>
fromPickle(std::istream &pickle);

//!  a pointer to the EnumerationStrategyBase
RDKIT_CHEMREACTIONS_EXPORT boost::shared_ptr<EnumerationStrategyBase>
fromPickle(const std::string &pickle);
}  // namespace EnumerationStrategyPickler
}  // namespace RDKit

#endif
