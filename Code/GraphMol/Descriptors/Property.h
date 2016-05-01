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
#include <boost/shared_ptr.hpp>

namespace RDKit {
namespace Descriptors {
  
// Computes an additive property from a molecule
struct PropertyFxn {
  virtual ~PropertyFxn() {}
  
  //! Compute the value for the molecule
  virtual double compute(const RDKit::ROMol &mol) const = 0;
  
  //! Is the property additive with respect to molecular fragments?
  virtual bool isAdditive() const { return false; }
};

//! Defines an additive property used to filter molecules
/* Known properties include
   MW - Molecular weight
   TPSA - polar surface area
   ALOGP - atom based log P
   NumRotors - the number of rotors (including side chain bond)
   MissingStereo - the number of atoms with unassigned stereo
   LipinskiHBA
   LipinskiHBD
   HBA
   HBD
   NumHeteroAtoms
   NumAmideBonds
   FractionCSP3
   NumRings
   NumAromaticRings
   NumAliphaticRings
   NumSaturatedRings
   NumHeteroCycles
   NumSaturatedHeteroCycles
   NumAromaticHeteroCycles
   NumAliphaticCarboCycles
   NumAromaticCarboCycles
   NumSaturatedCarboCycles
   NumAlpiphaticHeteroCycles
   NumSpiroAtoms
   NumBridgeheadAtoms
   AMW
   LabuteASA

   user defined types can be added as well with an appropriate
    computing functor.

   Example:
   \verbatim
     Property mw(Property::MW);
     Property logp(Property::ALOGP);

     // User based function
     struct MyPropertyFunction : public PropertyFxn {
       double compute(const ROMol &mol) const {
         return rdcast<double>(mol.getNumAtoms());
       }
       bool isAdditive() const { return true; }
     }

     // Make the a num atom property
     boost::shared_ptr<PropertyFxn> atomfxn(new MyPropertyFxn);
     Property numAtoms("NumAtoms", atomfxn);

     double alogp = logp.computeProperty(mol);
     double nAtoms = numAtoms.computeProperty(mol);
     
    \endverbatim

*/

class Property {
public:  
  enum PropType{ MW=0,
                 TPSA=1,
                 ALOGP=2,
                 NumRotors=3,
                 MissingStereo=4,
                 LipinskiHBA=5,
                 LipinskiHBD=6,
                 HBA=7,
                 HBD=8,
                 NumHeteroAtoms=9,
                 NumAmideBonds=10,
                 FractionCSP3=11,
                 NumRings=12,
                 NumAromaticRings=13,
                 NumAliphaticRings=14,
                 NumSaturatedRings=15,
                 NumHeteroCycles=16,
                 NumSaturatedHeteroCycles=17,
                 NumAromaticHeteroCycles=18,
                 NumAliphaticCarboCycles=19,
                 NumAromaticCarboCycles=21,
                 NumSaturatedCarboCycles=22,
                 NumAlpiphaticHeteroCycles=23,
                 NumSpiroAtoms=24,
                 NumBridgeheadAtoms=25,
                 AMW=26,
                 LabuteASA=27,
                 User=10000 // requires computing function
  };
  
private:  
  PropType    m_proptype;
  std::string m_propname;
  boost::shared_ptr<PropertyFxn> m_propfxn;

public:
  
 Property(PropType prop) :
  m_proptype(prop), m_propname(getNameForProp(prop)) {
  }

 Property(const std::string &name, boost::shared_ptr<PropertyFxn> fxn) :
  m_proptype(Property::User), m_propname(name), m_propfxn(fxn) {
  }
  
  // Descripts and versions for known property types
  static const char * getNameForProp(PropType prop);
  static std::string  getVersionForProp(PropType prop);

  const std::string &getName() const { return m_propname; }

  //! compute the property of a molecule
  double computeProperty(const ROMol &mol) const;

  //! Returns true if the property is additive with respect to fragments
  bool   isAdditive() const;
};

//! Defines a PropertyFilter used to filter molecules
//!  Using a minimum and a maximum value.
/* Known properties include
   MW - Molecular weight
   TPSA - polar surface area
   ALOGP - atom based log P
   NumRotors - the number of rotors (including side chain bond)
   MissingStereo - the number of atoms with unassigned stereo

   user defined types can be added as well with an appropriate
    computing functor.

   Example:
   \verbatim
     // Molecular weight filter 0 ... 500.
     PropertyFilter mw(PropertyFilter::MW, 0., 500.);
     
     // log p Filter -5 ... 5
     PropertyFilter logp(PropertyFilter::ALOGP, -5.0, 5.0);

     // User based function
     struct MyPropertyFunction : public PropertyFxn {
       double compute(const ROMol &mol) const {
         return rdcast<double>(mol.getNumAtoms());
       }
     }

     // add num atom property 0. ... 50.
     boost::shared_ptr<PropertyFxn> atomfxn(new MyPropertyFxn);
     PropertyFilter natoms("NumAtoms", 0., 50., atomfxn);

     if (natoms.accepts(mol)) {
      // molecule is ok
     } else {
      // molecules is bad
     }

    // These filters can also be used in FilterCatalogs
    //  example here...
    \endverbatim
*/
class PropertyFilter : public Property {
  double      m_min;
  double      m_max;
public:
  
PropertyFilter(PropType prop, double minValue, double maxValue) :  
  Property(prop), 
      m_min(minValue), m_max(maxValue) {
    PRECONDITION(prop != User, "User properties require a computing function");
  }
  
PropertyFilter(const std::string &propname,
               boost::shared_ptr<PropertyFxn> fxn,
               double minValue, double maxValue) :
  Property(propname, fxn),
      m_min(minValue), m_max(maxValue) {}
  
  //! returns the minimum acceptable value for the property
  double getMin() const { return m_min; }
  //! returns the maximum acceptable value for the property
  double getMax() const { return m_max; }

  //! returns true if the molecule passes the property filter
  bool accepts(const ROMol &mol) {
    double res = computeProperty(mol);
    return (res >= m_min && res <= m_max);
  }
};

}
}
#endif
