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
#include "Property.h"
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace RDKit {
namespace Descriptors {

namespace {
void _registerDescriptors() {
  REGISTER_DESCRIPTOR(exactmw, calcExactMW);
  REGISTER_DESCRIPTOR(amw, calcAMW);
  REGISTER_DESCRIPTOR(lipinskiHBA, calcLipinskiHBA);
  REGISTER_DESCRIPTOR(lipinskiHBD, calcLipinskiHBD);
  REGISTER_DESCRIPTOR(NumRotatableBonds, calcNumRotatableBonds);
  REGISTER_DESCRIPTOR(NumHBD, calcNumHBD);
  REGISTER_DESCRIPTOR(NumHBA, calcNumHBA);
  REGISTER_DESCRIPTOR(NumHeavyAtoms, calcNumHeavyAtoms);
  REGISTER_DESCRIPTOR(NumAtoms, calcNumAtoms);
  REGISTER_DESCRIPTOR(NumHeteroatoms, calcNumHeteroatoms);
  REGISTER_DESCRIPTOR(NumAmideBonds, calcNumAmideBonds);
  REGISTER_DESCRIPTOR(FractionCSP3, calcFractionCSP3);
  REGISTER_DESCRIPTOR(NumRings, calcNumRings);
  REGISTER_DESCRIPTOR(NumAromaticRings, calcNumAromaticRings);
  REGISTER_DESCRIPTOR(NumAliphaticRings, calcNumAliphaticRings);
  REGISTER_DESCRIPTOR(NumSaturatedRings, calcNumSaturatedRings);
  REGISTER_DESCRIPTOR(NumHeterocycles, calcNumHeterocycles);
  REGISTER_DESCRIPTOR(NumAromaticHeterocycles, calcNumAromaticHeterocycles);
  REGISTER_DESCRIPTOR(NumSaturatedHeterocycles, calcNumSaturatedHeterocycles);
  REGISTER_DESCRIPTOR(NumAliphaticHeterocycles, calcNumAliphaticHeterocycles);
  REGISTER_DESCRIPTOR(NumSpiroAtoms, calcNumSpiroAtoms);
  REGISTER_DESCRIPTOR(NumBridgeheadAtoms, calcNumBridgeheadAtoms);
  REGISTER_DESCRIPTOR(NumAtomStereoCenters, numAtomStereoCenters);
  REGISTER_DESCRIPTOR(NumUnspecifiedAtomStereoCenters,
                      numUnspecifiedAtomStereoCenters);
  REGISTER_DESCRIPTOR(labuteASA, calcLabuteASA);
  REGISTER_DESCRIPTOR(tpsa, calcTPSA);
  REGISTER_DESCRIPTOR(CrippenClogP, calcClogP);
  REGISTER_DESCRIPTOR(CrippenMR, calcMR);
  REGISTER_DESCRIPTOR(chi0v, calcChi0v);
  REGISTER_DESCRIPTOR(chi1v, calcChi1v);
  REGISTER_DESCRIPTOR(chi2v, calcChi3v);
  REGISTER_DESCRIPTOR(chi3v, calcChi3v);
  REGISTER_DESCRIPTOR(chi4v, calcChi4v);
  REGISTER_DESCRIPTOR(chi0n, calcChi0n);
  REGISTER_DESCRIPTOR(chi1n, calcChi1n);
  REGISTER_DESCRIPTOR(chi2n, calcChi3n);
  REGISTER_DESCRIPTOR(chi3n, calcChi3n);
  REGISTER_DESCRIPTOR(chi4n, calcChi4n);
  REGISTER_DESCRIPTOR(hallKierAlpha, calcHallKierAlpha);
  REGISTER_DESCRIPTOR(kappa1, calcKappa1);
  REGISTER_DESCRIPTOR(kappa2, calcKappa2);
  REGISTER_DESCRIPTOR(kappa3, calcKappa3);
  REGISTER_DESCRIPTOR(Phi, calcPhi);
};
}  // namespace

void registerDescriptors() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static std::once_flag once;
  std::call_once(once, _registerDescriptors);
#else
  static bool initialized = false;
  if (!initialized) {
    _registerDescriptors();
    initialized = true;
  }
#endif
}

std::vector<boost::shared_ptr<PropertyFunctor>> Properties::registry;
int Properties::registerProperty(PropertyFunctor *prop) {
  for (size_t i = 0; i < Properties::registry.size(); ++i) {
    if (registry[i]->getName() == prop->getName()) {
      Properties::registry[i] = boost::shared_ptr<PropertyFunctor>(prop);
      return i;
    }
  }
  // XXX Add mutex?
  Properties::registry.emplace_back(prop);
  return Properties::registry.size();
}

std::vector<std::string> Properties::getAvailableProperties() {
  registerDescriptors();
  std::vector<std::string> names;
  for (auto prop : Properties::registry) {
    names.push_back(prop->getName());
  }
  return names;
}

boost::shared_ptr<PropertyFunctor> Properties::getProperty(
    const std::string &name) {
  registerDescriptors();
  for (auto prop : Properties::registry) {
    if (prop.get() && prop->getName() == name) {
      return prop;
    }
  }
  throw KeyErrorException(name);
}

Properties::Properties() : m_properties() {
  registerDescriptors();
  for (auto prop : Properties::registry) {
    m_properties.push_back(prop);
  }
}

Properties::Properties(const std::vector<std::string> &propNames) {
  registerDescriptors();
  for (const auto &name : propNames) {
    m_properties.push_back(Properties::getProperty(name));
  }
}

std::vector<std::string> Properties::getPropertyNames() const {
  std::vector<std::string> names;
  for (auto prop : m_properties) {
    names.push_back(prop->getName());
  }
  return names;
}

std::vector<double> Properties::computeProperties(const RDKit::ROMol &mol,
                                                  bool annotate) const {
  std::vector<double> res;
  res.reserve(m_properties.size());
  for (auto prop : m_properties) {
    res.push_back((*prop)(mol));
    if (annotate) {
      mol.setProp<double>(prop->getName(), (*prop)(mol));
    }
  }
  return res;
}

void Properties::annotateProperties(RDKit::ROMol &mol) const {
  for (auto prop : m_properties) {
    mol.setProp<double>(prop->getName(), (*prop)(mol));
  }
}

PROP_RANGE_QUERY *makePropertyRangeQuery(const std::string &name, double min,
                                         double max) {
  auto *filter = new PROP_RANGE_QUERY(min, max);
  filter->setDataFunc(Properties::getProperty(name)->d_dataFunc);
  return filter;
}
}  // namespace Descriptors
}  // namespace RDKit
