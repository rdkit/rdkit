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
#include <RDGeneral/export.h>
#ifndef RDKIT_PROPERTIES_H
#define RDKIT_PROPERTIES_H

#include <GraphMol/RDKitBase.h>
#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <Query/Query.h>
#include <RDGeneral/Exceptions.h>

namespace RDKit {
namespace Descriptors {
struct RDKIT_DESCRIPTORS_EXPORT PropertyFunctor {
  // Registry of property functions
  //  See REGISTER_DESCRIPTOR
  std::string propName;
  std::string propVersion;
  double (*d_dataFunc)(const ROMol &);

  PropertyFunctor(const std::string &name, const std::string &version,
                  double (*func)(const ROMol &) = nullptr)
      : propName(name), propVersion(version), d_dataFunc(func) {}
  virtual ~PropertyFunctor(){};

  //! Compute the value of the property
  virtual double operator()(const RDKit::ROMol &) const = 0;

  //! Return the name of the property
  const std::string getName() const { return propName; }
  //! Return the properties version
  const std::string getVersion() const { return propVersion; }
};

//! Holds a collection of properties for computation purposes
class RDKIT_DESCRIPTORS_EXPORT Properties {
 protected:
  std::vector<boost::shared_ptr<PropertyFunctor>> m_properties;

 public:
  Properties();
  Properties(const std::vector<std::string> &propNames);

  std::vector<std::string> getPropertyNames() const;
  std::vector<double> computeProperties(const RDKit::ROMol &mol,
                                        bool annotate = false) const;
  void annotateProperties(RDKit::ROMol &mol) const;

  //! Register a property function - takes ownership
  static int registerProperty(PropertyFunctor *ptr);
  static boost::shared_ptr<PropertyFunctor> getProperty(
      const std::string &name);
  static std::vector<std::string> getAvailableProperties();
  static std::vector<boost::shared_ptr<PropertyFunctor>> registry;
};

typedef Queries::Query<bool, const ROMol &, true> PROP_BOOL_QUERY;
typedef Queries::AndQuery<int, const ROMol &, true> PROP_AND_QUERY;
typedef Queries::OrQuery<int, const ROMol &, true> PROP_OR_QUERY;
typedef Queries::XOrQuery<int, const ROMol &, true> PROP_XOR_QUERY;

typedef Queries::EqualityQuery<double, const ROMol &, true> PROP_EQUALS_QUERY;

typedef Queries::GreaterQuery<double, const ROMol &, true> PROP_GREATER_QUERY;

typedef Queries::GreaterEqualQuery<double, const ROMol &, true>
    PROP_GREATEREQUAL_QUERY;

typedef Queries::LessQuery<double, const ROMol &, true> PROP_LESS_QUERY;

typedef Queries::LessEqualQuery<double, const ROMol &, true>
    PROP_LESSEQUAL_QUERY;

typedef Queries::RangeQuery<double, const ROMol &, true> PROP_RANGE_QUERY;

template <class T>
T *makePropertyQuery(const std::string &name, double what) {
  T *t = new T(what);
  t->setDataFunc(Properties::getProperty(name)->d_dataFunc);
  return t;
}

RDKIT_DESCRIPTORS_EXPORT PROP_RANGE_QUERY *makePropertyRangeQuery(
    const std::string &name, double min, double max);

}  // namespace Descriptors
}  // namespace RDKit
#endif
