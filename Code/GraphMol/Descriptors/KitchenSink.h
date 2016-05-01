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
#ifndef RDKIT_KITCHENSINK_H
#define RDKIT_KITCHENSINK_H
#include "Property.h"
#include <vector>

namespace RDKit {
namespace Descriptors {

  const Property AllProperties[] = {
    Property::MW,
    Property::TPSA,
    Property::ALOGP,
    Property::NumRotors,
    Property::MissingStereo,
    Property::LipinskiHBA,
    Property::LipinskiHBD,
    Property::HBA,
    Property::HBD,
    Property::NumHeteroAtoms,
    Property::NumAmideBonds,
    Property::FractionCSP3,
    Property::NumRings,
    Property::NumAromaticRings,
    Property::NumAliphaticRings,
    Property::NumSaturatedRings,
    Property::NumHeteroCycles,
    Property::NumSaturatedHeteroCycles,
    Property::NumAromaticHeteroCycles,
    Property::NumAliphaticCarboCycles,
    Property::NumAromaticCarboCycles,
    Property::NumSaturatedCarboCycles,
    Property::NumAlpiphaticHeteroCycles,
    Property::NumSpiroAtoms,
    Property::NumBridgeheadAtoms,
    Property::AMW,
    Property::LabuteASA
  };

  const size_t NumProperties = sizeof(AllProperties)/sizeof(AllProperties[0]);
  
//! Holds a collection of properties for computation purposes
class Properties {
protected:
  std::vector<Property> m_properties;
  
public:
  Properties() : m_properties() {}
  Properties(const std::vector<Property> &props);
  
  std::vector<std::string> getPropertyNames() const;
  std::vector<double>      getProperties(const RDKit::ROMol &mol) const;
};

//! Computes all RDKit Properties
class KitchenSink : public Properties {
public:
  KitchenSink();
};

}
}
#endif
