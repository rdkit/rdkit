//  Copyright (c) 2017-2026, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <RDBoost/boost_shared_ptr.h>
#include <RDBoost/Wrap_nb.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/PatternFactory.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
std::string serializeSubstructLibrary(const SubstructLibrary &self) {
  return self.Serialize();
}

using GeneralizedSubstruct::ExtendedQueryMol;

const char *MolHolderBaseDoc =
    R"DOC(Base class for holding molecules used in the Substructure Library.
Instantiations of this class are passed into the SubstructureLibrary.
The API is quite simple:
  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
  GetMol(idx) -> return the molecule at index idx)DOC";

const char *MolHolderDoc =
    R"DOC(Holds raw in-memory molecules
  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
  GetMol(idx,sanitize=True) -> return the molecule at index idx)DOC";

const char *CachedMolHolderDoc =
    R"DOC(Holds molecules in their binary representation.
This allows more molecules to be held in memory at a time
  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

  AddBinary(data) -> adds a picked molecule molecule to the molecule holder, returns index of molecule
                     The data is stored as-is, no checking is done for validity.
  GetMol(idx) -> return the molecule at index idx)DOC";

const char *CachedSmilesMolHolderDoc =
    R"DOC(Holds molecules as smiles string
This allows more molecules to be held in memory at a time
  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

  AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                       The smiles is stored as-is, no checking is done for validity.
  GetMol(idx) -> return the molecule at index idx)DOC";

const char *CachedTrustedSmilesMolHolderDoc =
    R"DOC(Holds molecules as trusted smiles string
This allows more molecules to be held in memory at a time and avoids RDKit sanitization
overhead.
See: http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html
  AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

  AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                       The smiles is stored as-is, no checking is done for validity.
  GetMol(idx,s) -> return the molecule at index idx,
              note, only light sanitization is done here, for instance
              the molecules RingInfo is not initialized)DOC";

const char *PatternHolderDoc =
    R"DOC(Holds fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.)DOC";

const char *TautomerPatternHolderDoc =
    R"DOC(Holds tautomeric fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
These fingerprints are designed to be used with TautomerQueries.)DOC";

const char *KeyHolderDoc =
    R"DOC(Holds keys to return external references to the molecules in the molholder.
By default use the _Name property but can be overridden to be any property)DOC";

const char *SubstructLibraryDoc =
    R"DOC(SubstructLibrary: This provides a simple API for substructure searching large datasets
The SubstructLibrary takes full advantage of available threads during the search operation.
Basic operation is simple

>>> from __future__ import print_function
>>> import os
>>> from rdkit import Chem, RDConfig
>>> from rdkit.Chem import rdSubstructLibrary
>>> library = rdSubstructLibrary.SubstructLibrary()
>>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir,
...                               'NCI', 'first_200.props.sdf')):
...   idx = library.AddMol(mol)
>>> core = Chem.MolFromSmarts('CCCCOC')
>>> indices = library.GetMatches(core)
>>> len(indices)
11

Substructure matching options can be sent into GetMatches:

>>> indices = library.GetMatches(core, useChirality=False)
>>> len(indices)
11

Controlling the number of threads or the maximum number of matches returned:
is also available (the default is to run on all cores)

>>> indices = library.GetMatches(core, numThreads=2, maxResults=10)
>>> len(indices)
10

Working on larger datasets:

Molecules are fairly large objects and will limit the number that can be kept in memory.
To assist this we supply three other molecule holders:
  CachedMolHolder - stores molecules as their pickled representation

  CachedSmilesMolHolder - stores molecules internally as smiles strings

  CachedTrustedSmilesMolHolder = excepts (and stores) molecules as trusted smiles strings

Using Pattern fingerprints as a pre-filter:
Pattern fingerprints provide an easy way to indicate whether the substructure search should be
be done at all.  This is particularly useful with the Binary and Smiles based molecule holders
as they have an expensive molecule creation step in addition to the substructure searching step

>>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.CachedSmilesMolHolder(),
...                                               rdSubstructLibrary.PatternHolder())
>>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir,
...                               'NCI', 'first_200.props.sdf')):
...   idx = library.AddMol(mol)
>>> indices = library.GetMatches(core)
>>> len(indices)
11

This (obviously) takes longer to initialize.  However, both the molecule and pattern
holders can be populated with raw data, a simple example is below:

>>> import csv
>>> molholder = rdSubstructLibrary.CachedSmilesMolHolder()
>>> pattern_holder = rdSubstructLibrary.PatternHolder()
>>> with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.tpsa.csv')) as inf:
...   for i, row in enumerate(csv.reader(inf)):
...     if i:
...       idx = molholder.AddSmiles(row[0])
...       idx2 = pattern_holder.AddFingerprint(
...           pattern_holder.MakeFingerprint(Chem.MolFromSmiles(row[0])))
...       assert idx==idx2
>>> library = rdSubstructLibrary.SubstructLibrary(molholder,pattern_holder)
>>> indices = library.GetMatches(core)
>>> len(indices)
11

Finally, the KeyFromPropHolder can be used to use external keys such as
compound names.  By default the holder uses the '_Name' property but can
be changed to any property.

>>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.MolHolder(), rdSubstructLibrary.KeyFromPropHolder())
>>> m = Chem.MolFromSmiles('CCC')
>>> m.SetProp('_Name', 'Z11234')
>>> idx = library.AddMol(m)
>>> indices = library.GetMatches(m)
>>> list(library.GetKeyHolder().GetKeys(indices))
['Z11234']
)DOC";

nb::bytes SubstructLibrary_Serialize(const SubstructLibrary &self) {
  std::string res = self.Serialize();
  return nb::bytes(res.c_str(), res.size());
}

nb::tuple getSearchOrderHelper(const SubstructLibrary &self) {
  nb::list res;
  for (const auto v : self.getSearchOrder()) {
    res.append(v);
  }
  return nb::tuple(res);
}

void setSearchOrderHelper(SubstructLibrary &self, nb::handle seq) {
  nb::object seq_obj = nb::borrow(seq);
  std::unique_ptr<std::vector<unsigned int>> sorder =
      pythonObjectToVect<unsigned int>(seq_obj);
  if (sorder) {
    self.setSearchOrder(*sorder);
  } else {
    self.getSearchOrder().clear();
  }
}

void addPatternsHelper1(SubstructLibrary &self, int numThreads) {
  NOGIL gil;
  addPatterns(self, numThreads);
}

void addPatternsHelper2(SubstructLibrary &self,
                        boost::shared_ptr<FPHolderBase> patterns,
                        int numThreads) {
  NOGIL gil;
  addPatterns(self, patterns, numThreads);
}

// Template helper functions for the various Query types with GIL release
template <class Query>
std::vector<unsigned int> getMatchesGIL(const SubstructLibrary &self,
                                        const Query &query,
                                        bool recursionPossible,
                                        bool useChirality,
                                        bool useQueryQueryMatches,
                                        int numThreads, int maxResults) {
  NOGIL h;
  return self.getMatches(query, recursionPossible, useChirality,
                         useQueryQueryMatches, numThreads, maxResults);
}

template <class Query>
std::vector<unsigned int> getMatchesRangeGIL(
    const SubstructLibrary &self, const Query &query, unsigned int startIdx,
    unsigned int endIdx, bool recursionPossible, bool useChirality,
    bool useQueryQueryMatches, int numThreads, int maxResults) {
  NOGIL h;
  return self.getMatches(query, startIdx, endIdx, recursionPossible,
                         useChirality, useQueryQueryMatches, numThreads,
                         maxResults);
}

template <class Query>
std::vector<unsigned int> getMatchesParamsGIL(
    const SubstructLibrary &self, const Query &query,
    const SubstructMatchParameters &params, int numThreads, int maxResults) {
  NOGIL h;
  return self.getMatches(query, 0, self.size(), params, numThreads, maxResults);
}

template <class Query>
std::vector<unsigned int> getMatchesRangeParamsGIL(
    const SubstructLibrary &self, const Query &query, unsigned int startIdx,
    unsigned int endIdx, const SubstructMatchParameters &params, int numThreads,
    int maxResults) {
  NOGIL h;
  return self.getMatches(query, startIdx, endIdx, params, numThreads,
                         maxResults);
}

template <class Query>
unsigned int countMatchesGIL(const SubstructLibrary &self, const Query &query,
                             bool recursionPossible, bool useChirality,
                             bool useQueryQueryMatches, int numThreads) {
  NOGIL h;
  return self.countMatches(query, 0, self.size(), recursionPossible,
                           useChirality, useQueryQueryMatches, numThreads);
}

template <class Query>
unsigned int countMatchesRangeGIL(const SubstructLibrary &self,
                                  const Query &query, unsigned int startIdx,
                                  unsigned int endIdx, bool recursionPossible,
                                  bool useChirality, bool useQueryQueryMatches,
                                  int numThreads) {
  NOGIL h;
  return self.countMatches(query, startIdx, endIdx, recursionPossible,
                           useChirality, useQueryQueryMatches, numThreads);
}

template <class Query>
unsigned int countMatchesParamsGIL(const SubstructLibrary &self,
                                   const Query &query,
                                   const SubstructMatchParameters &params,
                                   int numThreads) {
  NOGIL h;
  return self.countMatches(query, 0, self.size(), params, numThreads);
}

template <class Query>
unsigned int countMatchesRangeParamsGIL(const SubstructLibrary &self,
                                        const Query &query,
                                        unsigned int startIdx,
                                        unsigned int endIdx,
                                        const SubstructMatchParameters &params,
                                        int numThreads) {
  NOGIL h;
  return self.countMatches(query, startIdx, endIdx, params, numThreads);
}

template <class Query>
bool hasMatchGIL(const SubstructLibrary &self, const Query &query,
                 bool recursionPossible, bool useChirality,
                 bool useQueryQueryMatches, int numThreads) {
  NOGIL h;
  return self.hasMatch(query, 0, self.size(), recursionPossible, useChirality,
                       useQueryQueryMatches, numThreads);
}

template <class Query>
bool hasMatchRangeGIL(const SubstructLibrary &self, const Query &query,
                      unsigned int startIdx, unsigned int endIdx,
                      bool recursionPossible, bool useChirality,
                      bool useQueryQueryMatches, int numThreads) {
  NOGIL h;
  return self.hasMatch(query, startIdx, endIdx, recursionPossible, useChirality,
                       useQueryQueryMatches, numThreads);
}

template <class Query>
bool hasMatchParamsGIL(const SubstructLibrary &self, const Query &query,
                       const SubstructMatchParameters &params, int numThreads) {
  NOGIL h;
  return self.hasMatch(query, 0, self.size(), params, numThreads);
}

template <class Query>
bool hasMatchRangeParamsGIL(const SubstructLibrary &self, const Query &query,
                            unsigned int startIdx, unsigned int endIdx,
                            const SubstructMatchParameters &params,
                            int numThreads) {
  NOGIL h;
  return self.hasMatch(query, startIdx, endIdx, params, numThreads);
}

const char *GetMatchesDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - numThreads: number of threads to use, -1 means all threads
  - maxResults: maximum number of results to return)DOC";

const char *GetMatchesRangeDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - startIdx:   index to search from
  - endIdx:     index (non-inclusize) to search to
  - numThreads: number of threads to use, -1 means all threads
  - maxResults: maximum number of results to return)DOC";

const char *CountMatchesDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - numThreads: number of threads to use, -1 means all threads)DOC";

const char *CountMatchesRangeDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - startIdx:   index to search from
  - endIdx:     index (non-inclusize) to search to
  - numThreads: number of threads to use, -1 means all threads)DOC";

const char *HasMatchDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - numThreads: number of threads to use, -1 means all threads)DOC";

const char *HasMatchRangeDoc =
    R"DOC(Get the matches for the query.

 Arguments:
  - query:      substructure query
  - startIdx:   index to search from
  - endIdx:     index (non-inclusize) to search to
  - numThreads: number of threads to use, -1 means all threads)DOC";

// clang-format off
#define LARGE_DEF(_tname_)                                                       \
  .def("GetMatches",                                                             \
       getMatchesGIL<_tname_>,                                                   \
       "query"_a,                                                                \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       "maxResults"_a = 1000,                                                    \
       GetMatchesDoc)                                                             \
  .def("GetMatches",                                                             \
       getMatchesRangeGIL<_tname_>,                                              \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       "maxResults"_a = 1000,                                                    \
       GetMatchesRangeDoc)                                                       \
  .def("GetMatches",                                                             \
       getMatchesParamsGIL<_tname_>,                                             \
       "query"_a,                                                                \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       "maxResults"_a = 1000,                                                    \
       GetMatchesDoc)                                                            \
  .def("GetMatches",                                                             \
       getMatchesRangeParamsGIL<_tname_>,                                        \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       "maxResults"_a = 1000,                                                    \
       GetMatchesRangeDoc)                                                       \
  .def("CountMatches",                                                           \
       countMatchesGIL<_tname_>,                                                 \
       "query"_a,                                                                \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       CountMatchesDoc)                                                          \
  .def("CountMatches",                                                           \
       countMatchesRangeGIL<_tname_>,                                            \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       CountMatchesRangeDoc)                                                     \
  .def("CountMatches",                                                           \
       countMatchesParamsGIL<_tname_>,                                           \
       "query"_a,                                                                \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       CountMatchesDoc)                                                          \
  .def("CountMatches",                                                           \
       countMatchesRangeParamsGIL<_tname_>,                                      \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       CountMatchesRangeDoc)                                                     \
  .def("HasMatch",                                                               \
       hasMatchGIL<_tname_>,                                                     \
       "query"_a,                                                                \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       HasMatchDoc)                                                              \
  .def("HasMatch",                                                               \
       hasMatchRangeGIL<_tname_>,                                                \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "recursionPossible"_a = true,                                             \
       "useChirality"_a = true,                                                  \
       "useQueryQueryMatches"_a = false,                                         \
       "numThreads"_a = -1,                                                      \
       HasMatchRangeDoc)                                                         \
  .def("HasMatch",                                                               \
       hasMatchParamsGIL<_tname_>,                                               \
       "query"_a,                                                                \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       HasMatchDoc)                                                              \
  .def("HasMatch",                                                               \
       hasMatchRangeParamsGIL<_tname_>,                                          \
       "query"_a,                                                                \
       "startIdx"_a,                                                             \
       "endIdx"_a,                                                               \
       "parameters"_a,                                                           \
       "numThreads"_a = -1,                                                      \
       HasMatchRangeDoc)
// clang-format on

void wrap_substructlibrary(nb::module_ &m) {
  nb::class_<MolHolderBase>(m, "MolHolderBase", MolHolderBaseDoc)
      .def("__len__", &MolHolderBase::size)
      .def("AddMol", &MolHolderBase::addMol, "m"_a,
           "Adds molecule to the molecule holder")
      .def("GetMol", &MolHolderBase::getMol, "arg1"_a,
           R"DOC(Returns a particular molecule in the molecule holder

  ARGUMENTS:
    - idx: which molecule to return

    - sanitize: if sanitize is False, return the internal molecule state [default True]

  NOTE: molecule indices start at 0)DOC");

  nb::class_<MolHolder, MolHolderBase>(m, "MolHolder", MolHolderDoc)
      .def(nb::init<>());

  nb::class_<CachedMolHolder, MolHolderBase>(m, "CachedMolHolder",
                                             CachedMolHolderDoc)
      .def(nb::init<>())
      .def(
          "AddBinary",
          [](CachedMolHolder &self, nb::bytes pickle) {
            return self.addBinary(std::string(
                static_cast<const char *>(pickle.data()), pickle.size()));
          },
          "pickle"_a,
          "Add a binary pickle to the molecule holder, no checking is done "
          "on the input data");

  nb::class_<CachedSmilesMolHolder, MolHolderBase>(m, "CachedSmilesMolHolder",
                                                   CachedSmilesMolHolderDoc)
      .def(nb::init<>())
      .def("AddSmiles", &CachedSmilesMolHolder::addSmiles, "smiles"_a,
           "Add a trusted smiles string to the molecule holder, no checking "
           "is done on the input data");

  nb::class_<CachedTrustedSmilesMolHolder, MolHolderBase>(
      m, "CachedTrustedSmilesMolHolder", CachedTrustedSmilesMolHolderDoc)
      .def(nb::init<>())
      .def("AddSmiles", &CachedTrustedSmilesMolHolder::addSmiles, "smiles"_a,
           "Add a trusted smiles string to the molecule holder, no checking "
           "is done on the input data");

  nb::class_<FPHolderBase>(m, "FPHolderBase", "")
      .def("__len__", &FPHolderBase::size)
      .def("AddMol", &FPHolderBase::addMol, "m"_a,
           "Adds a molecule to the fingerprint database, returns the index "
           "of the new pattern")
      .def("AddFingerprint",
           nb::overload_cast<const ExplicitBitVect &>(
               &FPHolderBase::addFingerprint),
           "v"_a,
           "Adds a raw bit vector to the fingerprint database, returns the "
           "index of the supplied pattern")
      .def("GetFingerprint", &FPHolderBase::getFingerprint,
           nb::rv_policy::reference_internal, "idx"_a,
           "Return the bit vector at the specified index")
      .def("PassesFilter", &FPHolderBase::passesFilter, "idx"_a, "query"_a,
           "Returns True if the specified index passes the filter supplied "
           "by the query bit vector")
      .def("MakeFingerprint", &FPHolderBase::makeFingerprint, "mol"_a,
           nb::rv_policy::take_ownership,
           "Compute the query bits for the holder");

  nb::class_<PatternHolder, FPHolderBase>(m, "PatternHolder", PatternHolderDoc)
      .def(nb::init<>())
      .def(nb::init<unsigned int>(), "numBits"_a);

  nb::class_<TautomerPatternHolder, FPHolderBase>(m, "TautomerPatternHolder",
                                                  TautomerPatternHolderDoc)
      .def(nb::init<>())
      .def(nb::init<unsigned int>(), "numBits"_a);

  nb::class_<KeyHolderBase>(m, "KeyHolderBase", "")
      .def("__len__", &KeyHolderBase::size)
      .def("AddMol", &KeyHolderBase::addMol, "m"_a,
           "Adds a molecule to the fingerprint database, returns the index "
           "of the new pattern")
      .def("AddKey", &KeyHolderBase::addKey, "arg1"_a,
           "Add a key to the key holder, must be manually synced")
      .def("GetKey", &KeyHolderBase::getKey, nb::rv_policy::copy, "arg1"_a,
           "Return the key at the specified index")
      .def("GetKeys", &KeyHolderBase::getKeys, "indices"_a,
           R"DOC(Returns the keys for the given indices as return by GetMatches

  ARGUMENTS:
    - indices: The indices of the keys)DOC");

  nb::class_<KeyFromPropHolder, KeyHolderBase>(m, "KeyFromPropHolder",
                                               KeyHolderDoc)
      .def(nb::init<>())
      .def(nb::init<const std::string &>(), "propname"_a)
      .def("GetPropName",
           nb::overload_cast<>(&KeyFromPropHolder::getPropName, nb::const_),
           nb::rv_policy::copy, "Return the key for the given molecule index");

  nb::class_<SubstructLibrary>(m, "SubstructLibrary", SubstructLibraryDoc)
      .def(nb::new_([]() { return new SubstructLibrary(); }))
      .def(nb::new_([](nb::bytes b) {
        return new SubstructLibrary(
            std::string(static_cast<const char *>(b.data()), b.size()));
      }))
      .def(nb::new_([](nb::handle mols_h) {
             MolHolderBase *mols = nb::cast<MolHolderBase *>(mols_h);
             return new SubstructLibrary(
                 nb::detail::shared_from_python<MolHolderBase>(mols, mols_h));
           }),
           "molecules"_a)
      .def(nb::new_([](nb::handle mols_h, nb::handle second) {
             // The second arg can be fingerprints (FPHolderBase) or keys
             // (KeyHolderBase) or None.
             MolHolderBase *mols = nb::cast<MolHolderBase *>(mols_h);
             boost::shared_ptr<MolHolderBase> molPtr =
                 nb::detail::shared_from_python<MolHolderBase>(mols, mols_h);
             if (second.is_none()) {
               return new SubstructLibrary(molPtr);
             } else if (nb::isinstance<FPHolderBase>(second)) {
               FPHolderBase *fps = nb::cast<FPHolderBase *>(second);
               return new SubstructLibrary(
                   molPtr,
                   nb::detail::shared_from_python<FPHolderBase>(fps, second));
             } else {
               KeyHolderBase *keys = nb::cast<KeyHolderBase *>(second);
               return new SubstructLibrary(
                   molPtr,
                   nb::detail::shared_from_python<KeyHolderBase>(keys, second));
             }
           }),
           "molecules"_a, "fingerprints_or_keys"_a.none())
      .def(nb::new_([](nb::handle mols_h, nb::handle fps_h, nb::handle keys_h) {
             MolHolderBase *mols = nb::cast<MolHolderBase *>(mols_h);
             boost::shared_ptr<MolHolderBase> molPtr =
                 nb::detail::shared_from_python<MolHolderBase>(mols, mols_h);
             boost::shared_ptr<FPHolderBase> fpHolder;
             if (!fps_h.is_none()) {
               FPHolderBase *fps = nb::cast<FPHolderBase *>(fps_h);
               fpHolder =
                   nb::detail::shared_from_python<FPHolderBase>(fps, fps_h);
             }
             boost::shared_ptr<KeyHolderBase> keyHolder;
             if (!keys_h.is_none()) {
               KeyHolderBase *keys = nb::cast<KeyHolderBase *>(keys_h);
               keyHolder =
                   nb::detail::shared_from_python<KeyHolderBase>(keys, keys_h);
             }
             return new SubstructLibrary(molPtr, fpHolder, keyHolder);
           }),
           "molecules"_a, "fingerprints"_a.none(), "keys"_a.none())

      .def("GetMolHolder",
           [](SubstructLibrary &self) { return self.getMolHolder(); })
      .def("GetFpHolder",
           [](SubstructLibrary &self) { return self.getFpHolder(); })
      .def("GetKeyHolder",
           [](SubstructLibrary &self) { return self.getKeyHolder(); })

      .def("AddMol", &SubstructLibrary::addMol, "mol"_a,
           "Adds a molecule to the substruct library")

      // clang-format off
      LARGE_DEF(ROMol)
      LARGE_DEF(TautomerQuery)
      LARGE_DEF(MolBundle)
      LARGE_DEF(ExtendedQueryMol)
      // clang-format on

      .def("GetMol", &SubstructLibrary::getMol, "idx"_a,
           R"DOC(Returns a particular molecule in the molecule holder

  ARGUMENTS:
    - idx: which molecule to return

  NOTE: molecule indices start at 0)DOC")

      .def("SetSearchOrder", setSearchOrderHelper, "seq"_a.none(),
           R"DOC(Sets the search order for the library

  ARGUMENTS:
    - order: sequence of molecule indices

  NOTE: molecule indices start at 0)DOC")
      .def("GetSearchOrder", getSearchOrderHelper,
           R"DOC(Returns the search order for the library

  NOTE: molecule indices start at 0)DOC")

      .def("__len__", &SubstructLibrary::size)

      .def(
          "ToStream",
          [](const SubstructLibrary &self, nb::object &fileobj) {
            std::ostringstream oss;
            self.toStream(oss);
            fileobj.attr("write")(oss.str());
          },
          "stream"_a,
          R"DOC(Serialize a substructure library to a python text stream.
The stream can be a file in text mode or an io.StringIO type object

  ARGUMENTS:
    - stream: a text or text stream like object

  >>> from rdkit.Chem import rdSubstructLibrary
  >>> import io
  >>> lib = rdSubstructLibrary.SubstructLibrary()
  >>> stream = io.StringIO()
  >>> lib.ToStream(stream)

   or
  >>> with open('rdkit.sslib', 'w') as stream:
  ...  lib.ToStream(stream))DOC")

      .def(
          "InitFromStream",
          [](SubstructLibrary &self, nb::object &fileobj) {
            nb::object data = fileobj.attr("read")();
            std::string s;
            if (nb::isinstance<nb::bytes>(data)) {
              auto b = nb::cast<nb::bytes>(data);
              s = std::string(static_cast<const char *>(b.data()), b.size());
            } else {
              s = nb::cast<std::string>(data);
            }
            std::istringstream iss(s);
            self.initFromStream(iss);
          },
          "stream"_a,
          R"DOC(Deserialize a substructure library from a python bytes stream.
Python doesn't allow seeking operations inside a unicode or string stream anymore
so this requires opening a file in binary mode or using an io.ByteIO type object

  ARGUMENTS:
    - stream: a binary stream like object

  SubstructLibrary.Serialize already writes a binary stream

  >>> from rdkit.Chem import rdSubstructLibrary
  >>> import io
  >>> lib = rdSubstructLibrary.SubstructLibrary()
  >>> stream = io.BytesIO( lib.Serialize() )
  >>> lib.InitFromStream(stream)

   remember to write to text and read from a binary stream
  >>> with open('rdkit.sslib', 'w') as f: lib.ToStream(f)
  >>> with open('rdkit.sslib', 'rb') as f: lib.InitFromStream(f))DOC")

      .def("Serialize", SubstructLibrary_Serialize)
      .def("__getstate__",
           getObjectState<SubstructLibrary, serializeSubstructLibrary>)
      .def("__setstate__", setObjectState<SubstructLibrary>);

  m.def("SubstructLibraryCanSerialize", SubstructLibraryCanSerialize,
        "Returns True if the SubstructLibrary is serializable "
        "(requires boost serialization");

  m.def("AddPatterns", addPatternsHelper1, "sslib"_a, "numThreads"_a = 1,
        "Add pattern fingerprints to the given library, use "
        "numThreads=-1 to use all available cores");

  m.def("AddPatterns", addPatternsHelper2, "sslib"_a, "patterns"_a,
        "numThreads"_a = 1,
        "Add pattern fingerprints to the given library, use numThreads=-1 to "
        "use all available cores");
}

}  // namespace RDKit
