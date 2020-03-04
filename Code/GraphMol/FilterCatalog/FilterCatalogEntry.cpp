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

#include "Filters.h"
#include "FilterCatalogEntry.h"
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit {

const std::string DESCRIPTION = "description";
std::string FilterCatalogEntry::getDescription() const {
  std::string desc;
  d_props.getValIfPresent(DESCRIPTION, desc);
  return desc;
}

void FilterCatalogEntry::setDescription(const std::string &desc) {
  d_props.setVal(DESCRIPTION, desc);
}

void FilterCatalogEntry::toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#endif
}

std::string FilterCatalogEntry::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void FilterCatalogEntry::initFromStream(std::istream &ss) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#endif
}

void FilterCatalogEntry::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

FilterCatalogEntry *MakeFilterCatalogEntry(const FilterData_t &data,
                                           unsigned int num_props,
                                           const FilterProperty_t *props) {
  const int debugParse = 0;
  const bool mergeHs = true;
  ROMOL_SPTR pattern(SmartsToMol(data.smarts, debugParse, mergeHs));
  if (!pattern.get()) {
    return nullptr;
  }

  // The filter has the concept of the maximum number of times the pattern is
  // allowed
  //   without triggering
  // The matcher has the concept of the minimum pattern count before the pattern
  //  is triggered.
  // Hence pattern.minCount == filter.maxCount + 1
  const unsigned int minCount = data.max ? data.max + 1 : 1;
  FilterCatalogEntry *entry = new FilterCatalogEntry(
      data.name, boost::shared_ptr<FilterMatcherBase>(
                     new SmartsMatcher(data.name, pattern, minCount)));

  if (!entry->isValid()) {
    delete entry;
    return nullptr;
  }

  // add the props
  for (unsigned int i = 0; i < num_props; i++) {
    std::string val(props[i].value);
    entry->getProps().setVal<std::string>(std::string(props[i].key), val);
  }

  return entry;
}
}  // end namespace RDKit
