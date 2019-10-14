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
#ifndef __RD_FILTER_CATALOG_H__
#define __RD_FILTER_CATALOG_H__

#include <RDGeneral/types.h>  // For Dict
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Catalogs/CatalogEntry.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

#include "FilterMatchers.h"

namespace RDKit {
typedef std::map<std::string, std::string> STRING_PROPS;

class RDKIT_FILTERCATALOG_EXPORT FilterCatalogEntry
    : public RDCatalog::CatalogEntry {
 private:
  boost::shared_ptr<FilterMatcherBase> d_matcher;
  Dict d_props;

 public:
  FilterCatalogEntry() : d_matcher(), d_props() {}

  FilterCatalogEntry(const std::string &name, const FilterMatcherBase &matcher)
      : RDCatalog::CatalogEntry(), d_matcher(matcher.copy()) {
    setDescription(name);
  }

  FilterCatalogEntry(const std::string &name,
                     boost::shared_ptr<FilterMatcherBase> matcher)
      : RDCatalog::CatalogEntry(), d_matcher(matcher) {
    setDescription(name);
  }

  FilterCatalogEntry(const FilterCatalogEntry &rhs)
      : RDCatalog::CatalogEntry(rhs),
        d_matcher(rhs.d_matcher),
        d_props(rhs.d_props) {}

  virtual ~FilterCatalogEntry() {}

  //------------------------------------
  //! Returns true if the Filters stored in this catalog entry are valid

  bool isValid() const { return d_matcher.get() && d_matcher->isValid(); }

  //------------------------------------
  //! Returns the description of the catalog entry
  std::string getDescription() const;

  //------------------------------------
  //! Sets the description of the catalog entry
  /*
    \param const std::string & description
   */
  void setDescription(const std::string &description);

  //! \name Properties
  //@{

  //! returns a list with the names of our \c properties
  STR_VECT getPropList() const { return d_props.keys(); }

  //! sets a \c property value
  /*!
     \param key the name under which the \c property should be stored.
         If a \c property is already stored under this name, it will be
         replaced.
     \param val the value to be stored
   */
  template <typename T>
  void setProp(const char *key, T val) {
    std::string what(key);
    setProp(what, val);
  }
  //! \overload
  template <typename T>
  void setProp(const std::string &key, T val) {
    d_props.setVal(key, val);
  }

  //! allows retrieval of a particular property value
  /*!

     \param key the name under which the \c property should be stored.
         If a \c property is already stored under this name, it will be
         replaced.
     \param res a reference to the storage location for the value.

     <b>Notes:</b>
       - if no \c property with name \c key exists, a KeyErrorException will be
     thrown.
       - the \c boost::lexical_cast machinery is used to attempt type
     conversions.
         If this fails, a \c boost::bad_lexical_cast exception will be thrown.

  */
  template <typename T>
  void getProp(const char *key, T &res) const {
    d_props.getVal(key, res);
  }
  //! \overload
  template <typename T>
  void getProp(const std::string &key, T &res) const {
    d_props.getVal(key, res);
  }
  //! \overload
  template <typename T>
  T getProp(const char *key) const {
    return d_props.getVal<T>(key);
  }
  //! \overload
  template <typename T>
  T getProp(const std::string &key) const {
    return d_props.getVal<T>(key);
  }

  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  template <typename T>
  bool getPropIfPresent(const char *key, T &res) const {
    return d_props.getValIfPresent(key, res);
  }
  //! \overload
  template <typename T>
  bool getPropIfPresent(const std::string &key, T &res) const {
    return d_props.getValIfPresent(key, res);
  }

  //! returns whether or not we have a \c property with name \c key
  bool hasProp(const char *key) const { return d_props.hasVal(key); }
  //! \overload
  bool hasProp(const std::string &key) const {
    return d_props.hasVal(key);
    // return hasProp(key.c_str());
  }

  //! clears the value of a \c property
  void clearProp(const char *key) {
    std::string what(key);
    clearProp(what);
  };
  //! \overload
  void clearProp(const std::string &key) { d_props.clearVal(key); };

  // -------------------------------------------
  //!  Properties usually contain the reference and source
  //!  for the catalog entry.

  Dict &getProps() { return d_props; }
  const Dict &getProps() const { return d_props; }
  void setProps(const Dict &props) { d_props = props; }

  //------------------------------------
  //! Returns the matching filters for this catalog entry
  /*
    \param mol The molecule to match against
    \param std::vector<FilterMatch> the resulting FilterMatches
   */
  bool getFilterMatches(const ROMol &mol,
                        std::vector<FilterMatch> &matchVect) const {
    return this->isValid() && d_matcher->getMatches(mol, matchVect);
  }

  //------------------------------------
  //! Returns true if the filters in this catalog entry match the molecule
  /*
    \param mol The molecule to match against
   */

  bool hasFilterMatch(const ROMol &mol) const {
    return this->isValid() && d_matcher->hasMatch(mol);
  }

  //! serializes (pickles) to a stream
  virtual void toStream(std::ostream &ss) const;
  //! returns a string with a serialized (pickled) representation
  virtual std::string Serialize() const;
  //! initializes from a stream pickle
  virtual void initFromStream(std::istream &ss);
  //! initializes from a string pickle
  virtual void initFromString(const std::string &text);

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    RDUNUSED_PARAM(version);
    registerFilterMatcherTypes(ar);

    ar &d_matcher;
    // we only save string based props here...
    STR_VECT keys = d_props.keys();
    std::vector<std::string> string_props;
    for (size_t i = 0; i < keys.size(); ++i) {
      std::string val;
      try {
        if (d_props.getValIfPresent<std::string>(keys[i], val)) {
          string_props.push_back(keys[i]);
          string_props.push_back(val);
        }
      } catch (const boost::bad_any_cast &) {
        // pass, can't serialize
        // warning, this changes properties types, see Dict.cpp
      }
    }
    ar &string_props;
  }

  template <class Archive>
  void load(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    registerFilterMatcherTypes(ar);

    ar &d_matcher;
    std::vector<std::string> string_props;
    ar &string_props;
    d_props.reset();

    for (size_t i = 0; i < string_props.size() / 2; ++i) {
      d_props.setVal(string_props[i * 2], string_props[i * 2 + 1]);
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();
#endif
};
}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::FilterCatalogEntry, 1);
#endif

#endif  //__RD_FILTER_CATALOG_H__
