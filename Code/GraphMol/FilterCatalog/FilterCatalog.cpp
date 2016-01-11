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

#include "FilterCatalog.h"
#include "Filters.h"
#include "FilterMatchers.h"

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#endif

#include <climits>  // to get CHAR_BIT # bits in a char

namespace RDKit {
bool FilterCatalogParams::addCatalog(FilterCatalogs catalog) {
  bool addedCatalog = false;
  // meh, not the best here... perhaps make num_catalogs?
  for (unsigned int i = 0; i < sizeof(catalog) * CHAR_BIT; ++i) {
    if ((catalog & (1u << i)) & FilterCatalogParams::ALL) {
      FilterCatalogs cat = static_cast<FilterCatalogs>(catalog & (1u << i));
      if (GetNumEntries(cat)) {
        d_catalogs.push_back(cat);
        addedCatalog = true;
      }
    }
  }

  return addedCatalog;
}

void FilterCatalogParams::fillCatalog(FilterCatalog &catalog) {
  for (size_t i = 0; i < getCatalogs().size(); ++i) {
    const FilterCatalogs catalogToAdd = getCatalogs()[i];

    const unsigned int entries = GetNumEntries(catalogToAdd);
    const unsigned int propEntries = GetNumPropertyEntries(catalogToAdd);
    // XXX Fix Me -> these should probably be shared to save memory
    const FilterProperty_t *props = GetFilterProperties(catalogToAdd);

    for (unsigned int i = 0; i < entries; ++i) {
      const FilterData_t &data = GetFilterData(catalogToAdd)[i];
      FilterCatalogEntry *entry =
          MakeFilterCatalogEntry(data, propEntries, props);
      PRECONDITION(entry, "Bad Entry data");
      if (entry) catalog.addEntry(entry);  // catalog owns entry
    }
  }
}

void FilterCatalogParams::toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#endif
}

std::string FilterCatalogParams::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void FilterCatalogParams::initFromStream(std::istream &ss) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#endif
}

void FilterCatalogParams::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

FilterCatalog::~FilterCatalog() {}

void FilterCatalog::Clear() { d_entries.clear(); }

FilterCatalog::FilterCatalog(const std::string &binStr)
    : FCatalog(), d_entries() {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  std::stringstream ss(binStr);
  boost::archive::text_iarchive ar(ss);
  ar &d_entries;
#endif
}

std::string FilterCatalog::Serialize() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else

  std::stringstream ss;
  boost::archive::text_oarchive ar(ss);
  ar &d_entries;
  return ss.str();
#endif
}

unsigned int FilterCatalog::addEntry(entryType_t *entry, bool updateFPLength) {
  RDUNUSED_PARAM(updateFPLength);
  return addEntry(boost::shared_ptr<entryType_t>(entry));
}

unsigned int FilterCatalog::addEntry(SENTRY entry, bool updateFPLength) {
  RDUNUSED_PARAM(updateFPLength);
  d_entries.push_back(entry);
  return static_cast<unsigned int>(d_entries.size() - 1);
}

const FilterCatalog::entryType_t *FilterCatalog::getEntryWithIdx(
    unsigned int idx) const {
  if (idx < d_entries.size()) return d_entries[idx].get();
  return 0;
}

FilterCatalog::CONST_SENTRY FilterCatalog::getEntry(unsigned int idx) const {
  PRECONDITION(idx < d_entries.size(), "Index out of bounds");
  if (idx < d_entries.size()) return d_entries[idx];
  return CONST_SENTRY();
}

bool FilterCatalog::removeEntry(unsigned int idx) {
  if (idx < d_entries.size()) {
    d_entries.erase(d_entries.begin() + (idx));
    return true;
  }
  return false;
}

bool FilterCatalog::removeEntry(const FilterCatalog::CONST_SENTRY &entry) {
  std::vector<SENTRY>::iterator it =
      std::find(d_entries.begin(), d_entries.end(), entry);
  if (it != d_entries.end()) {
    d_entries.erase(it);
    return true;
  }
  return false;
}

unsigned int FilterCatalog::getIdxForEntry(const entryType_t *entry) const {
  for (size_t i = 0; i < d_entries.size(); ++i) {
    if (d_entries[i].get() == entry) return i;
  }
  return UINT_MAX;
}

unsigned int FilterCatalog::getIdxForEntry(const CONST_SENTRY &entry) const {
  for (size_t i = 0; i < d_entries.size(); ++i) {
    if (d_entries[i] == entry) return i;
  }
  return UINT_MAX;
}

void FilterCatalog::setCatalogParams(paramType_t *params) {
  Clear();
  FCatalog::setCatalogParams(params);
  params->fillCatalog(*this);
}

bool FilterCatalog::hasMatch(const ROMol &mol) const {
  return getFirstMatch(mol) != 0;
}

FilterCatalog::CONST_SENTRY FilterCatalog::getFirstMatch(
    const ROMol &mol) const {
  for (size_t i = 0; i < d_entries.size(); ++i) {
    if (d_entries[i]->hasFilterMatch(mol)) return d_entries[i];
  }
  return CONST_SENTRY();
}

const std::vector<FilterCatalog::CONST_SENTRY> FilterCatalog::getMatches(
    const ROMol &mol) const {
  std::vector<CONST_SENTRY> result;
  for (size_t i = 0; i < d_entries.size(); ++i) {
    if (d_entries[i]->hasFilterMatch(mol)) result.push_back(d_entries[i]);
  }
  return result;
}

bool FilterCatalogCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
}
}
