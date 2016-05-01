//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
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

#include "KitchenSink.h"
#include "Property.h"

namespace RDKit {
namespace Descriptors {


Properties::Properties(const std::vector<Property> &props) :
    m_properties(props) {}

std::vector<std::string> Properties::getPropertyNames() const {
  std::vector<std::string> res;
  res.reserve(m_properties.size());
  for (size_t i=0; i<m_properties.size(); ++i) {
    res.push_back(m_properties[i].getName());
  }
  return res;
}

std::vector<double> Properties::getProperties(const RDKit::ROMol &mol) const {
  std::vector<double> res;
  res.reserve(m_properties.size());
  for (size_t i=0; i<m_properties.size(); ++i) {
    res.push_back(m_properties[i].computeProperty(mol));
  }
  return res;
}

KitchenSink::KitchenSink() :
    Properties() {
  m_properties = std::vector<Property>(
      AllProperties,
      AllProperties + sizeof(AllProperties)/sizeof(AllProperties[0]));
}
}
};
