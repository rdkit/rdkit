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
#ifndef __RD_FILTER_MATCHER_BASE_H__
#define __RD_FILTER_MATCHER_BASE_H__
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif  // RDK_USE_BOOST_SERIALIZATION

namespace RDKit {

class FilterMatcherBase;  // Forward declaration

//! Holds the atomPairs matched by the underlying matcher
struct RDKIT_FILTERCATALOG_EXPORT FilterMatch {
  boost::shared_ptr<FilterMatcherBase> filterMatch;
  MatchVectType atomPairs;

  FilterMatch() : filterMatch(), atomPairs() {}
  FilterMatch(boost::shared_ptr<FilterMatcherBase> filter,
              MatchVectType atomPairs)
      : filterMatch(filter), atomPairs(atomPairs) {}

  FilterMatch(const FilterMatch &rhs)
      : filterMatch(rhs.filterMatch), atomPairs(rhs.atomPairs) {}

  bool operator==(const FilterMatch &rhs) const {
    return (filterMatch.get() == rhs.filterMatch.get() &&
            atomPairs == rhs.atomPairs);
  }

  bool operator!=(const FilterMatch &rhs) const {
    return !(filterMatch.get() == rhs.filterMatch.get() &&
             atomPairs == rhs.atomPairs);
  }
};

RDKIT_FILTERCATALOG_EXPORT extern const char *DEFAULT_FILTERMATCHERBASE_NAME;
class RDKIT_FILTERCATALOG_EXPORT FilterMatcherBase
    : public boost::enable_shared_from_this<FilterMatcherBase> {
  //------------------------------------
  //! Virtual API for filter matching
  std::string d_filterName;

 public:
  FilterMatcherBase(const std::string &name = DEFAULT_FILTERMATCHERBASE_NAME)
      : boost::enable_shared_from_this<FilterMatcherBase>(),
        d_filterName(name) {}

  FilterMatcherBase(const FilterMatcherBase &rhs)
      : boost::enable_shared_from_this<FilterMatcherBase>(),
        d_filterName(rhs.d_filterName) {}

  virtual ~FilterMatcherBase() {}

  virtual bool isValid() const = 0;

  virtual std::string getName() const { return d_filterName; }
  //------------------------------------
  //! getMatches
  /*!
    Match the filter against a molecule

    \param mol readonly const molecule being searched
    \param matches  output vector of atom index matches found in the molecule
  */

  virtual bool getMatches(const ROMol &mol,
                          std::vector<FilterMatch> &matchVect) const = 0;

  //------------------------------------
  //! hasMatches
  /*!
    Does the given molecule contain this filter pattern

    \param mol readonly const molecule being searched
  */

  virtual bool hasMatch(const ROMol &mol) const = 0;

  //------------------------------------
  //! Clone - deprecated
  //  Clones the current FilterMatcherBase into one that
  //   can be passed around safely.
  virtual boost::shared_ptr<FilterMatcherBase> Clone() const {
    BOOST_LOG(rdWarningLog)
        << "FilterMatcherBase::Clone is deprecated, use copy instead"
        << std::endl;
    return copy();
  }

  //------------------------------------
  //! copy
  //  copies the current FilterMatcherBase into one that
  //   can be passed around safely.
  virtual boost::shared_ptr<FilterMatcherBase> copy() const = 0;

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &d_filterName;
  }
#endif
};

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_SERIALIZATION_ASSUME_ABSTRACT(FilterMatcherBase)
#endif
}  // namespace RDKit
#endif
