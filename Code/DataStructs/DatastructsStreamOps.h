//  Copyright (c) 2019, Novartis Institutes for BioMedical Research Inc.
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
#ifndef RDKIT_DATASTRUCTS_STREAMOPS
#define RDKIT_DATASTRUCTS_STREAMOPS
#include <RDGeneral/StreamOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <typeinfo>
#include <boost/any.hpp>

namespace RDKit {
class DataStructsExplicitBitVecPropHandler : public CustomPropHandler {
 public:
  const char *getPropName() const { return "ExplicitBVProp"; }
  bool canSerialize(const RDValue &value) const {
    return rdvalue_is<ExplicitBitVect>(value);
  }

  bool read(std::istream &ss, RDValue &value) const {
    std::string v;
    int version = 0;
    streamRead(ss, v, version);
    ExplicitBitVect bv(v);
    value = bv;
    return true;
  }

  bool write(std::ostream &ss, const RDValue &value) const {
    try {
      std::string output =
          rdvalue_cast<const ExplicitBitVect &>(value).toString();
      streamWrite(ss, output);
      return true;
    } catch (boost::bad_any_cast &) {
      return false;
    }
  }

  CustomPropHandler *clone() const {
    return new DataStructsExplicitBitVecPropHandler;
  }
};

}  // namespace RDKit
#endif
