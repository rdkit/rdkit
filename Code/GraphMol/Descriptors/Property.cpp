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

#ifdef RDK_THREADSAFE_SSS
#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread/once.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {
namespace Descriptors {

namespace {
const std::string CrippenClogPVersion = crippenVersion;
double calcClogP(const ROMol &mol) {
  double clogp,mr;
  calcCrippenDescriptors(mol, clogp, mr);
  return clogp;
}
const std::string CrippenMRVersion = crippenVersion;
double calcMR(const ROMol &mol) {
  double clogp,mr;
  calcCrippenDescriptors(mol, clogp, mr);
  return mr;
}
}

namespace {
void _registerDescriptors() {
  REGISTER_DESCRIPTOR(exactmw, calcExactMW);
  REGISTER_DESCRIPTOR(lipinskiHBA, calcLipinskiHBA);
  REGISTER_DESCRIPTOR(lipinskiHBD, calcLipinskiHBD);
  REGISTER_DESCRIPTOR(NumRotatableBonds, calcNumRotatableBonds);
  REGISTER_DESCRIPTOR(NumHBD, calcNumHBD);
  REGISTER_DESCRIPTOR(NumHBA, calcNumHBA);
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
  REGISTER_DESCRIPTOR(NumUnspecifiedAtomStereoCenters, numUnspecifiedAtomStereoCenters);
  REGISTER_DESCRIPTOR(labuteASA, calcLabuteASA);
  REGISTER_DESCRIPTOR(tpsa, calcTPSA);
  REGISTER_DESCRIPTOR(CrippenClogP, calcClogP);
  REGISTER_DESCRIPTOR(CrippenMR, calcMR);
};
}

void registerDescriptors() {
#ifdef RDK_THREADSAFE_SSS
  static boost::once_flag once = BOOST_ONCE_INIT;
  boost::call_once(&_registerDescriptors, once);
#else
  static bool initialized = false;
  if(!initialized) {
    _registerDescriptors();
    initalized=true;
  }
#endif
}

std::vector<boost::shared_ptr<PropertyFxn> > Properties::registry;
int Properties::registerProperty(PropertyFxn *prop) {
  for(size_t i=0; i<Properties::registry.size(); ++i) {
    if (registry[i]->getName() == prop->getName()) {
      Properties::registry[i] = boost::shared_ptr<PropertyFxn>(prop);
      return i;
    }
  }
  // XXX Add mutex?
  Properties::registry.push_back( boost::shared_ptr<PropertyFxn>(prop) );
  return Properties::registry.size();
}

std::vector<std::string> Properties::getAvailableProperties() {
  registerDescriptors();
  std::vector<std::string> names;
  BOOST_FOREACH(boost::shared_ptr<PropertyFxn> prop, Properties::registry) {
    names.push_back(prop->getName());
  }
  return names;
}

boost::shared_ptr<PropertyFxn> Properties::getProperty(const std::string &name) {
  registerDescriptors();
  BOOST_FOREACH(boost::shared_ptr<PropertyFxn> prop, Properties::registry) {
    if (prop.get() && prop->getName() == name) {
      return prop;
    }
  }
  throw KeyErrorException(name);
}

Properties::Properties() : m_properties() {
  registerDescriptors();
  BOOST_FOREACH(boost::shared_ptr<PropertyFxn> prop, Properties::registry) {
    m_properties.push_back(prop);
  }
}

Properties::Properties(const std::vector<std::string> &propNames) {
  registerDescriptors();
  BOOST_FOREACH(const std::string &name, propNames) {
    m_properties.push_back( Properties::getProperty(name) );
  }
}

std::vector<std::string> Properties::getPropertyNames() const {
  std::vector<std::string> names;
  BOOST_FOREACH(boost::shared_ptr<PropertyFxn> prop, m_properties) {
    names.push_back(prop->getName());
  }
  return names;
}

std::vector<double> Properties::computeProperties(const RDKit::ROMol &mol) const {
  std::vector<double> res;
  res.reserve(m_properties.size());
  BOOST_FOREACH(boost::shared_ptr<PropertyFxn> prop, m_properties) {
    res.push_back(prop->compute(mol));
  }
  return res;
}

}
}
