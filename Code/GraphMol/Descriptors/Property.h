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
#ifndef RDKIT_PROPERTIES_H
#define RDKIT_PROPERTIES_H

#include <GraphMol/RDKitBase.h>
#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace Descriptors {
  struct PropertyFxn {
    // Registry of property functions
    //  See REGISTER_DESCRIPTOR
    std::string propName;
    std::string propVersion;
    
   PropertyFxn(const std::string &name, const std::string &version) :
    propName(name), propVersion(version) {
    }
    virtual ~PropertyFxn() {};
    
    //! Compute the value of the property
    virtual double compute(const RDKit::ROMol &) const = 0;
    
    //! Return the name of the property
    const std::string getName() const { return propName; }
    //! Return the properties version
    const std::string getVersion() const { return propVersion; }
    
};

  
//! Holds a collection of properties for computation purposes
class Properties {
protected:
  std::vector<boost::shared_ptr<PropertyFxn> > m_properties;
  
public:
  Properties();
  Properties(const std::vector<std::string> &propNames);
  
  std::vector<std::string> getPropertyNames() const;
  std::vector<double>      computeProperties(const RDKit::ROMol &mol) const;

  //! Register a property function - takes ownership
  static int registerProperty(PropertyFxn *ptr);
  static boost::shared_ptr<PropertyFxn> getProperty(const std::string &name);
  static std::vector<std::string> getAvailableProperties();
  static std::vector<boost::shared_ptr<PropertyFxn> > registry;
    
  
};

//! Filters properties between a min and max
class PropertyFilter {

  boost::shared_ptr<PropertyFxn> m_property;
  double      m_min;
  double      m_max;
  
public:
  
  PropertyFilter(const std::string &propName, double minValue, double maxValue) :  
   m_property(Properties::getProperty(propName)), m_min(minValue), m_max(maxValue) {
  }
  
  PropertyFilter(boost::shared_ptr<PropertyFxn> fxn,
               double minValue, double maxValue) :
   m_property(fxn), m_min(minValue), m_max(maxValue) {}
  
  //! returns the minimum acceptable value for the property
  double getMin() const { return m_min; }
  //! returns the maximum acceptable value for the property
  double getMax() const { return m_max; }

  //! returns true if the molecule passes the property filter
  bool accepts(const ROMol &mol) {
    double res = m_property->compute(mol);
    return (res >= m_min && res <= m_max);
  }
};



}
}
#endif
