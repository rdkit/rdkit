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
#ifndef _RD_FILTER_CATALOG_PARAMS_
#define _RD_FILTER_CATALOG_PARAMS_

#include <Catalogs/Catalog.h>
#include <Catalogs/CatalogParams.h>
#include "FilterCatalogEntry.h"

namespace RDKit {
class FilterCatalog;
class RDKIT_FILTERCATALOG_EXPORT FilterCatalogParams
    : public RDCatalog::CatalogParams {
 public:
  enum FilterCatalogs {
    PAINS_A = (1u << 1),
    PAINS_B = (1u << 2),
    PAINS_C = (1u << 3),
    PAINS = PAINS_A | PAINS_B | PAINS_C,

    BRENK = (1u << 4),
    NIH = (1u << 5),
    ZINC = (1u << 6),

    ALL = PAINS | BRENK | NIH | ZINC
  };

  FilterCatalogParams() : RDCatalog::CatalogParams() {
    setTypeStr("Filter Catalog Parameters");
  }

  FilterCatalogParams(FilterCatalogs catalogs) : RDCatalog::CatalogParams() {
    setTypeStr("Filter Catalog Parameters");
    addCatalog(catalogs);
  }

  FilterCatalogParams(const FilterCatalogParams &other)
      : RDCatalog::CatalogParams(other), d_catalogs(other.d_catalogs) {}

  virtual ~FilterCatalogParams() {}

  //------------------------------------
  //! Adds an existing FilterCatalog specification to be used in the
  // FilterCatalog
  //
  /*!
    Specifies an existing filter catalog to be used.

      \param catalogs One of the enumerated known FilterCatalogs
    */
  virtual bool addCatalog(FilterCatalogs catalogs);

  //------------------------------------
  //! Returns the existing list of FilterCatalogs to be used.
  const std::vector<FilterCatalogs> &getCatalogs() const { return d_catalogs; }
  //! Fill a catalog with the appropriate entries
  virtual void fillCatalog(FilterCatalog &catalog) const;

  //! serializes (pickles) to a stream
  virtual void toStream(std::ostream &ss) const;
  //! returns a string with a serialized (pickled) representation
  virtual std::string Serialize() const;
  //! initializes from a stream pickle
  virtual void initFromStream(std::istream &ss);
  //! initializes from a string pickle
  virtual void initFromString(const std::string &text);

 private:
  std::vector<FilterCatalogs> d_catalogs;

#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &d_catalogs;
  }
#endif
};

typedef RDCatalog::Catalog<FilterCatalogEntry, FilterCatalogParams> FCatalog;
class RDKIT_FILTERCATALOG_EXPORT FilterCatalog : public FCatalog {
 public:
  // syntactic sugar for getMatch(es) return values.
  typedef boost::shared_ptr<FilterCatalogEntry> SENTRY;

  // If boost::python can support shared_ptr of const objects
  //  we can enable support for this feature
  typedef boost::shared_ptr<const entryType_t> CONST_SENTRY;

  FilterCatalog() : FCatalog(), d_entries() {}

  FilterCatalog(FilterCatalogParams::FilterCatalogs catalogs)
      : FCatalog(), d_entries() {
    paramType_t temp_params(catalogs);
    setCatalogParams(&temp_params);
  }

  FilterCatalog(const FilterCatalogParams &params) : FCatalog(), d_entries() {
    setCatalogParams(&params);
  }

  FilterCatalog(const FilterCatalog &rhs)
      : FCatalog(rhs), d_entries(rhs.d_entries) {}

  FilterCatalog(const std::string &binStr);

  ~FilterCatalog();

  virtual std::string Serialize() const;

  // Adds a new FilterCatalogEntry to the catalog
  /*!
    Adds a new FilterCatalogEntry to the catalog  The catalog
    owns the entry

    \param entry          The FilterCatalogEntry to add.
    \param updateFPLength unused in the FilterCatalog object.
  */

  virtual unsigned int addEntry(FilterCatalogEntry *entry,
                                bool updateFPLength = true);

  // Adds a new FilterCatalogEntry to the catalog
  /*!
    Adds a new FilterCatalogEntry to the catalog  The catalog
    owns the entry

    \param entry          The shared_ptr of the FilterCatalogEntry to add.
    \param updateFPLength unused in the FilterCatalog object.
  */

  virtual unsigned int addEntry(SENTRY entry, bool updateFPLength = true);

  // Removes a FilterCatalogEntry to the catalog by description
  /*!
    Removes a FilterCatalogEntry from the catalog.

    \param idx          The FilterCatalogEntry index for the entry to remove.
     n.b. removing an entry may change the indices of other entries.
          To safely remove entries, remove entries with the highest idx
           first.
  */
  bool removeEntry(unsigned int idx);
  bool removeEntry(CONST_SENTRY entry);

  //------------------------------------
  //! returns a particular FilterCatalogEntry in the Catalog
  //!  required by Catalog.h API
  virtual const FilterCatalogEntry *getEntryWithIdx(unsigned int idx) const;

  //------------------------------------
  //! returns a particular FilterCatalogEntry in the Catalog
  //!  memory safe version of getEntryWithIdx
  CONST_SENTRY getEntry(unsigned int idx) const;

  //------------------------------------
  //! returns the idx of the given entry, UINT_MAX if not found.

  unsigned int getIdxForEntry(const FilterCatalogEntry *entry) const;
  unsigned int getIdxForEntry(CONST_SENTRY entry) const;

  //------------------------------------
  //! returns the number of entries in the catalog
  virtual unsigned int getNumEntries() const {
    return static_cast<unsigned int>(d_entries.size());
  }

  //------------------------------------
  //! Reset the current catalog to match the specified FilterCatalogParameters
  /*
    \param params  The new FilterCatalogParams specifying the new state of the
    catalog
  */
  virtual void setCatalogParams(const FilterCatalogParams *params);

  //------------------------------------
  //! Returns true if the molecule matches any entry in the catalog
  /*
    \param mol  ROMol to match against the catalog
  */
  bool hasMatch(const ROMol &mol) const;

  //------------------------------------
  //! Returns the first match against the catalog
  /*
    \param mol  ROMol to match against the catalog
  */
  CONST_SENTRY getFirstMatch(const ROMol &mol) const;

  //-------------------------------------------
  //! Returns all entry matches to the molecule
  /*
    \param mol  ROMol to match against the catalog
  */
  const std::vector<CONST_SENTRY> getMatches(const ROMol &mol) const;

  //--------------------------------------------
  //! Returns all FilterMatches for the molecule
  /*
    \param mol  ROMol to match against the catalog
  */
  const std::vector<FilterMatch> getFilterMatches(const ROMol &mol) const;

 private:
  void Clear();
  std::vector<SENTRY> d_entries;
};

RDKIT_FILTERCATALOG_EXPORT bool FilterCatalogCanSerialize();

//! Run a filter catalog on a set of smiles strings
/*
  \param smiles vector of smiles strings to analyze
  \param nthreads specify the number of threads to use or specify 0 to use all processors
                         [default 1]
  \returns a vector of vectors.  For each input smiles string, returns
                   a vector of shared_ptr::FilterMatchEntry objects.
                   If a molecule matches no filters, the vector will be empty.
                   If a smiles can't be parsed, a 'no valid RDKit molecule' catalog entry is returned.

*/
RDKIT_FILTERCATALOG_EXPORT
std::vector<std::vector<boost::shared_ptr<const FilterCatalogEntry>>> RunFilterCatalog(
              const FilterCatalog &filterCatalog,
	      const std::vector<std::string> &smiles,
	      int numThreads=1);
}  // namespace RDKit

#endif
