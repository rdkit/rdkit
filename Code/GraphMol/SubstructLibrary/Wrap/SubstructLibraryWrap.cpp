//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
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
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/python_streambuf.h>

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/PatternFactory.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>

namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit {

using GeneralizedSubstruct::ExtendedQueryMol;

// Because we need to release the GIL before we launch a thread, we need to make
// a thin stub
//  for every function that does this.  This stub exists Because I couldn't
//  quite figure out
// how to make a release GIL call_guard to work with the complexity of the
// functions here.
//
// We could have made helper functions, but this was easier in the end as all we
// needed
//  to do was replace SubstructLibrary with SubstructLibraryWrap in the existing
//  boost::python wrap below.
class SubstructLibraryWrap {
 public:
  SubstructLibrary ss;

  SubstructLibraryWrap() : ss() {}
  SubstructLibraryWrap(boost::shared_ptr<MolHolderBase> molecules)
      : ss(molecules) {}
  SubstructLibraryWrap(boost::shared_ptr<MolHolderBase> molecules,
                       boost::shared_ptr<FPHolderBase> fingerprints)
      : ss(molecules, fingerprints) {}
  SubstructLibraryWrap(boost::shared_ptr<MolHolderBase> molecules,
                       boost::shared_ptr<KeyHolderBase> keys)
      : ss(molecules, keys) {}
  SubstructLibraryWrap(boost::shared_ptr<MolHolderBase> molecules,
                       boost::shared_ptr<FPHolderBase> fingerprints,
                       boost::shared_ptr<KeyHolderBase> keys)
      : ss(molecules, fingerprints, keys) {}
  SubstructLibraryWrap(const std::string &pickle) : ss(pickle) {}

  boost::shared_ptr<MolHolderBase> &getMolHolder() { return ss.getMolHolder(); }
  boost::shared_ptr<FPHolderBase> &getFpHolder() { return ss.getFpHolder(); }
  boost::shared_ptr<KeyHolderBase> &getKeyHolder() { return ss.getKeyHolder(); }
  unsigned int addMol(const ROMol &mol) { return ss.addMol(mol); }

  template <class Query>
  std::vector<unsigned int> getMatches(const Query &query,
                                       bool recursionPossible = true,
                                       bool useChirality = true,
                                       bool useQueryQueryMatches = false,
                                       int numThreads = -1,
                                       int maxResults = -1) const {
    NOGIL h;
    return ss.getMatches(query, recursionPossible, useChirality,
                         useQueryQueryMatches, numThreads, maxResults);
  }
  template <class Query>
  std::vector<unsigned int> getMatches(const Query &query,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const {
    NOGIL h;
    return ss.getMatches(query, 0, size(), params, numThreads, maxResults);
  }

  template <class Query>
  std::vector<unsigned int> getMatches(
      const Query &query, unsigned int startIdx, unsigned int endIdx,
      bool recursionPossible = true, bool useChirality = true,
      bool useQueryQueryMatches = false, int numThreads = -1,
      int maxResults = -1) const {
    NOGIL h;
    return ss.getMatches(query, startIdx, endIdx, recursionPossible,
                         useChirality, useQueryQueryMatches, numThreads,
                         maxResults);
  };

  template <class Query>
  std::vector<unsigned int> getMatches(const Query &query,
                                       unsigned int startIdx,
                                       unsigned int endIdx,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const {
    NOGIL h;
    return ss.getMatches(query, startIdx, endIdx, params, numThreads,
                         maxResults);
  }

  template <class Query>
  unsigned int countMatches(const Query &query, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1) const {
    NOGIL h;
    return ss.countMatches(query, 0, size(), recursionPossible, useChirality,
                           useQueryQueryMatches, numThreads);
  }

  template <class Query>
  unsigned int countMatches(const Query &query,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const {
    NOGIL h;
    return ss.countMatches(query, 0, size(), params, numThreads);
  }

  template <class Query>
  unsigned int countMatches(const Query &query, unsigned int startIdx,
                            unsigned int endIdx, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1) const {
    NOGIL h;
    return ss.countMatches(query, startIdx, endIdx, recursionPossible,
                           useChirality, useQueryQueryMatches, numThreads);
  };

  template <class Query>
  unsigned int countMatches(const Query &query, unsigned int startIdx,
                            unsigned int endIdx,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const {
    NOGIL h;
    return ss.countMatches(query, startIdx, endIdx, params, numThreads);
  }

  template <class Query>
  bool hasMatch(const Query &query, bool recursionPossible = true,
                bool useChirality = true, bool useQueryQueryMatches = false,
                int numThreads = -1) const {
    NOGIL h;
    return ss.hasMatch(query, 0, size(), recursionPossible, useChirality,
                       useQueryQueryMatches, numThreads);
  }

  template <class Query>
  bool hasMatch(const Query &query, const SubstructMatchParameters &params,
                int numThreads = -1) const {
    NOGIL h;
    return ss.hasMatch(query, 0, size(), params, numThreads);
  }
  template <class Query>
  bool hasMatch(const Query &query, unsigned int startIdx, unsigned int endIdx,
                bool recursionPossible = true, bool useChirality = true,
                bool useQueryQueryMatches = false, int numThreads = -1) const {
    NOGIL h;
    return ss.hasMatch(query, startIdx, endIdx, recursionPossible, useChirality,
                       useQueryQueryMatches, numThreads);
  };

  template <class Query>
  bool hasMatch(const Query &query, unsigned int startIdx, unsigned int endIdx,
                const SubstructMatchParameters &params,
                int numThreads = -1) const {
    NOGIL h;
    return ss.hasMatch(query, startIdx, endIdx, params, numThreads);
  }

  boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    return ss.getMol(idx);
  }
  unsigned int size() const { return ss.size(); }
};

const char *MolHolderBaseDoc =
    "Base class for holding molecules used in the Substructure Library.\n"
    "Instantiations of this class are passed into the SubstructureLibrary.\n"
    "The API is quite simple: \n"
    "  AddMol(mol) -> adds a molecule to the molecule holder, returns index of "
    "molecule\n"
    "  GetMol(idx) -> return the molecule at index idx\n";

const char *MolHolderDoc =
    "Holds raw in-memory molecules\n"
    "  AddMol(mol) -> adds a molecule to the molecule holder, returns index of "
    "molecule\n"
    "  GetMol(idx,sanitize=True) -> return the molecule at index idx\n";

const char *CachedMolHolderDoc =
    "Holds molecules in their binary representation.\n"
    "This allows more molecules to be held in memory at a time\n"
    "  AddMol(mol) -> adds a molecule to the molecule holder, returns index of "
    "molecule\n\n"
    "  AddBinary(data) -> adds a picked molecule molecule to the molecule "
    "holder, returns index of molecule\n"
    "                     The data is stored as-is, no checking is done for "
    "validity.\n"
    "  GetMol(idx) -> return the molecule at index idx\n";

const char *CachedSmilesMolHolderDoc =
    "Holds molecules as smiles string\n"
    "This allows more molecules to be held in memory at a time\n"
    "  AddMol(mol) -> adds a molecule to the molecule holder, returns index of "
    "molecule\n\n"
    "  AddSmiles(smiles) -> adds a smiles string to the molecule holder, "
    "returns index of molecule\n"
    "                       The smiles is stored as-is, no checking is done "
    "for validity.\n"
    "  GetMol(idx) -> return the molecule at index idx\n";

const char *CachedTrustedSmilesMolHolderDoc =
    "Holds molecules as trusted smiles string\n"
    "This allows more molecules to be held in memory at a time and avoids "
    "RDKit sanitization\n"
    "overhead.\n"
    "See: "
    "http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html\n"
    "  AddMol(mol) -> adds a molecule to the molecule holder, returns index of "
    "molecule\n\n"
    "  AddSmiles(smiles) -> adds a smiles string to the molecule holder, "
    "returns index of molecule\n"
    "                       The smiles is stored as-is, no checking is done "
    "for validity.\n"
    "  GetMol(idx,s) -> return the molecule at index idx, \n"
    "              note, only light sanitization is done here, for instance\n"
    "              the molecules RingInfo is not initialized\n";

const char *PatternHolderDoc =
    "Holds fingerprints with optional, user-defined number of bits (default: "
    "2048) used for filtering of molecules.";
const char *TautomerPatternHolderDoc =
    "Holds tautomeric fingerprints with optional, user-defined number of bits "
    "(default: "
    "2048) used for filtering of molecules.\n"
    "These fingerprints are designed to be used with TautomerQueries.";

const char *KeyHolderDoc =
    "Holds keys to return external references to the molecules in the "
    "molholder.\n"
    "By default use the _Name property but can be overridden to be any "
    "property";

const char *SubstructLibraryDoc =
    "SubstructLibrary: This provides a simple API for substructure searching "
    "large datasets\n"
    "The SubstructLibrary takes full advantage of available threads during the "
    "search operation.\n"
    "Basic operation is simple\n"
    "\n"
    ">>> from __future__ import print_function\n"
    ">>> import os\n"
    ">>> from rdkit import Chem, RDConfig\n"
    ">>> from rdkit.Chem import rdSubstructLibrary\n"
    ">>> library = rdSubstructLibrary.SubstructLibrary()\n"
    ">>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir, \n"
    "...                               'NCI', 'first_200.props.sdf')):\n"
    "...   idx = library.AddMol(mol)\n"
    ">>> core = Chem.MolFromSmarts('CCCCOC')\n"
    ">>> indices = library.GetMatches(core)\n"
    ">>> len(indices)\n"
    "11\n"
    "\n"
    "Substructure matching options can be sent into GetMatches:\n"
    "\n"
    ">>> indices = library.GetMatches(core, useChirality=False) \n"
    ">>> len(indices)\n"
    "11\n"
    "\n"
    "Controlling the number of threads or the maximum number of matches "
    "returned:\n"
    "is also available (the default is to run on all cores)\n"
    "\n"
    ">>> indices = library.GetMatches(core, numThreads=2, maxResults=10) \n"
    ">>> len(indices)\n"
    "10\n"
    "\n"
    "Working on larger datasets:\n"
    "\n"
    "Molecules are fairly large objects and will limit the number that can be "
    "kept in memory.\n"
    "To assist this we supply three other molecule holders:\n"
    "  CachedMolHolder - stores molecules as their pickled representation\n"
    "\n"
    "  CachedSmilesMolHolder - stores molecules internally as smiles strings\n"
    "\n"
    "  CachedTrustedSmilesMolHolder = excepts (and stores) molecules as "
    "trusted smiles strings\n"
    "\n"
    "Using Pattern fingerprints as a pre-filter:"
    "\n"
    "Pattern fingerprints provide an easy way to indicate whether the "
    "substructure search should be\n"
    "be done at all.  This is particularly useful with the Binary and Smiles "
    "based molecule holders\n"
    "as they have an expensive molecule creation step in addition to the "
    "substructure searching step\n "
    "\n"
    ">>> library = "
    "rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary."
    "CachedSmilesMolHolder(), \n"
    "...                                               "
    "rdSubstructLibrary.PatternHolder())\n"
    ">>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir, \n"
    "...                               'NCI', 'first_200.props.sdf')):\n"
    "...   idx = library.AddMol(mol)\n"
    ">>> indices = library.GetMatches(core)\n"
    ">>> len(indices)\n"
    "11\n"
    "\n"
    "This (obviously) takes longer to initialize.  However, both the molecule "
    "and pattern\n"
    "holders can be populated with raw data, a simple example is below:\n"
    "\n"
    ">>> import csv\n"
    ">>> molholder = rdSubstructLibrary.CachedSmilesMolHolder()\n"
    ">>> pattern_holder = rdSubstructLibrary.PatternHolder()\n"
    ">>> with open(os.path.join(RDConfig.RDDataDir, 'NCI', "
    "'first_200.tpsa.csv')) as inf:\n"
    "...   for i, row in enumerate(csv.reader(inf)):\n"
    "...     if i:\n"
    "...       idx = molholder.AddSmiles(row[0])\n"
    "...       idx2 = pattern_holder.AddFingerprint(\n"
    "...           "
    "pattern_holder.MakeFingerprint(Chem.MolFromSmiles(row[0])))\n"
    "...       assert idx==idx2\n"
    ">>> library = "
    "rdSubstructLibrary.SubstructLibrary(molholder,pattern_holder)\n"
    ">>> indices = library.GetMatches(core)\n"
    ">>> len(indices)\n"
    "11\n"
    "\n"
    "Finally, the KeyFromPropHolder can be used to use external keys such as\n"
    "compound names.  By default the holder uses the '_Name' property but can\n"
    "be changed to any property.\n"
    "\n"
    ">>> library = "
    "rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.MolHolder(), "
    "rdSubstructLibrary.KeyFromPropHolder())\n"
    ">>> m = Chem.MolFromSmiles('CCC')\n"
    ">>> m.SetProp('_Name', 'Z11234')\n"
    ">>> idx = library.AddMol(m)\n"
    ">>> indices = library.GetMatches(m)\n"
    ">>> list(library.GetKeyHolder().GetKeys(indices))\n"
    "['Z11234']\n"
    "";

python::object SubstructLibrary_Serialize(const SubstructLibraryWrap &cat) {
  std::string res = cat.ss.Serialize();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

struct substructlibrary_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const SubstructLibraryWrap &self) {
    std::string res;
    if (!SubstructLibraryCanSerialize()) {
      throw_runtime_error("Pickling of FilterCatalog instances is not enabled");
    }
    res = self.ss.Serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

void toStream(const SubstructLibraryWrap &cat, python::object &fileobj) {
  streambuf ss(fileobj, 't');
  streambuf::ostream ost(ss);
  cat.ss.toStream(ost);
}

void initFromStream(SubstructLibraryWrap &cat, python::object &fileobj) {
  streambuf ss(fileobj,
               'b');  // python StringIO can't seek, so need binary data
  streambuf::istream is(ss);
  cat.ss.initFromStream(is);
}

boost::shared_ptr<MolHolderBase> GetMolHolder(SubstructLibraryWrap &sslib) {
  // need to convert from a ref to a real shared_ptr
  return sslib.ss.getMolHolder();
}

boost::shared_ptr<FPHolderBase> GetFpHolder(SubstructLibraryWrap &sslib) {
  // need to convert from a ref to a real shared_ptr
  return sslib.ss.getFpHolder();
}

boost::shared_ptr<KeyHolderBase> GetKeyHolder(SubstructLibraryWrap &sslib) {
  // need to convert from a ref to a real shared_ptr
  return sslib.ss.getKeyHolder();
}

python::tuple getSearchOrderHelper(const SubstructLibraryWrap &sslib) {
  python::list res;
  for (const auto v : sslib.ss.getSearchOrder()) {
    res.append(v);
  }
  return python::tuple(res);
}
void setSearchOrderHelper(SubstructLibraryWrap &sslib,
                          const python::object &seq) {
  std::unique_ptr<std::vector<unsigned int>> sorder =
      pythonObjectToVect<unsigned int>(seq);
  if (sorder) {
    sslib.ss.setSearchOrder(*sorder);
  } else {
    sslib.ss.getSearchOrder().clear();
  }
}

void addPatternsHelper(SubstructLibraryWrap &sslib,
                       boost::shared_ptr<FPHolderBase> patterns,
                       int numThreads) {
  NOGIL gil;
  addPatterns(sslib.ss, patterns, numThreads);
}

void addPatternsHelper(SubstructLibraryWrap &sslib, int numThreads) {
  NOGIL gil;
  addPatterns(sslib.ss, numThreads);
}

#define LARGE_DEF(_tname_)                                                     \
  .def("GetMatches",                                                           \
       (std::vector<unsigned int>(SubstructLibraryWrap::*)(                    \
           const _tname_ &, bool, bool, bool, int, int) const) &               \
           SubstructLibraryWrap::getMatches,                                   \
       (python::arg("self"), python::arg("query"),                             \
        python::arg("recursionPossible") = true,                               \
        python::arg("useChirality") = true,                                    \
        python::arg("useQueryQueryMatches") = false,                           \
        python::arg("numThreads") = -1, python::arg("maxResults") = 1000),     \
       "Get the matches for the query.\n\n"                                    \
       " Arguments:\n"                                                         \
       "  - query:      substructure query\n"                                  \
       "  - numThreads: number of threads to use, -1 means all threads\n"      \
       "  - maxResults: maximum number of results to return")                  \
      .def(                                                                    \
          "GetMatches",                                                        \
          (std::vector<unsigned int>(SubstructLibraryWrap::*)(                 \
              const _tname_ &, unsigned int, unsigned int, bool, bool, bool,   \
              int, int) const) &                                               \
              SubstructLibraryWrap::getMatches,                                \
          (python::arg("self"), python::arg("query"), python::arg("startIdx"), \
           python::arg("endIdx"), python::arg("recursionPossible") = true,     \
           python::arg("useChirality") = true,                                 \
           python::arg("useQueryQueryMatches") = false,                        \
           python::arg("numThreads") = -1, python::arg("maxResults") = 1000),  \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - startIdx:   index to search from\n"                             \
          "  - endIdx:     index (non-inclusize) to search to\n"               \
          "  - numThreads: number of threads to use, -1 means all threads\n"   \
          "  - maxResults: maximum number of results to return")               \
      .def("CountMatches",                                                     \
           (unsigned int (SubstructLibraryWrap::*)(const _tname_ &, bool,      \
                                                   bool, bool, int) const) &   \
               SubstructLibraryWrap::countMatches,                             \
           (python::arg("self"), python::arg("query"),                         \
            python::arg("recursionPossible") = true,                           \
            python::arg("useChirality") = true,                                \
            python::arg("useQueryQueryMatches") = false,                       \
            python::arg("numThreads") = -1),                                   \
           "Get the matches for the query.\n\n"                                \
           " Arguments:\n"                                                     \
           "  - query:      substructure query\n"                              \
           "  - numThreads: number of threads to use, -1 means all threads\n") \
      .def(                                                                    \
          "CountMatches",                                                      \
          (unsigned int (SubstructLibraryWrap::*)(                             \
              const _tname_ &, unsigned int, unsigned int, bool, bool, bool,   \
              int) const) &                                                    \
              SubstructLibraryWrap::countMatches,                              \
          (python::arg("self"), python::arg("query"), python::arg("startIdx"), \
           python::arg("endIdx"), python::arg("recursionPossible") = true,     \
           python::arg("useChirality") = true,                                 \
           python::arg("useQueryQueryMatches") = false,                        \
           python::arg("numThreads") = -1),                                    \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - startIdx:   index to search from\n"                             \
          "  - endIdx:     index (non-inclusize) to search to\n"               \
          "  - numThreads: number of threads to use, -1 means all threads\n")  \
      .def("HasMatch",                                                         \
           (bool(SubstructLibraryWrap::*)(const _tname_ &, bool, bool, bool,   \
                                          int) const) &                        \
               SubstructLibraryWrap::hasMatch,                                 \
           (python::arg("self"), python::arg("query"),                         \
            python::arg("recursionPossible") = true,                           \
            python::arg("useChirality") = true,                                \
            python::arg("useQueryQueryMatches") = false,                       \
            python::arg("numThreads") = -1),                                   \
           "Get the matches for the query.\n\n"                                \
           " Arguments:\n"                                                     \
           "  - query:      substructure query\n"                              \
           "  - numThreads: number of threads to use, -1 means all threads\n") \
      .def(                                                                    \
          "HasMatch",                                                          \
          (bool(SubstructLibraryWrap::*)(const _tname_ &, unsigned int,        \
                                         unsigned int, bool, bool, bool, int)  \
               const) &                                                        \
              SubstructLibraryWrap::hasMatch,                                  \
          (python::arg("self"), python::arg("query"), python::arg("startIdx"), \
           python::arg("endIdx"), python::arg("recursionPossible") = true,     \
           python::arg("useChirality") = true,                                 \
           python::arg("useQueryQueryMatches") = false,                        \
           python::arg("numThreads") = -1),                                    \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - startIdx:   index to search from\n"                             \
          "  - endIdx:     index (non-inclusize) to search to\n"               \
          "  - numThreads: number of threads to use, -1 means all threads\n")  \
      .def("GetMatches",                                                       \
           (std::vector<unsigned int>(SubstructLibraryWrap::*)(                \
               const _tname_ &, const SubstructMatchParameters &, int, int)    \
                const) &                                                       \
               SubstructLibraryWrap::getMatches,                               \
           (python::arg("self"), python::arg("query"),                         \
            python::arg("parameters"), python::arg("numThreads") = -1,         \
            python::arg("maxResults") = 1000),                                 \
           "Get the matches for the query.\n\n"                                \
           " Arguments:\n"                                                     \
           "  - query:      substructure query\n"                              \
           "  - numThreads: number of threads to use, -1 means all threads\n"  \
           "  - maxResults: maximum number of results to return")              \
      .def(                                                                    \
          "GetMatches",                                                        \
          (std::vector<unsigned int>(SubstructLibraryWrap::*)(                 \
              const _tname_ &, unsigned int, unsigned int,                     \
              const SubstructMatchParameters &, int, int) const) &             \
              SubstructLibraryWrap::getMatches,                                \
          (python::arg("self"), python::arg("query"), python::arg("startIdx"), \
           python::arg("endIdx"), python::arg("parameters"),                   \
           python::arg("numThreads") = -1, python::arg("maxResults") = 1000),  \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - startIdx:   index to search from\n"                             \
          "  - endIdx:     index (non-inclusize) to search to\n"               \
          "  - numThreads: number of threads to use, -1 means all threads\n"   \
          "  - maxResults: maximum number of results to return")               \
      .def(                                                                    \
          "CountMatches",                                                      \
          (unsigned int (SubstructLibraryWrap::*)(                             \
              const _tname_ &, const SubstructMatchParameters &, int) const) & \
              SubstructLibraryWrap::countMatches,                              \
          (python::arg("self"), python::arg("query"),                          \
           python::arg("parameters"), python::arg("numThreads") = -1),         \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - numThreads: number of threads to use, -1 means all threads\n")  \
      .def("CountMatches",                                                     \
           (unsigned int (SubstructLibraryWrap::*)(                            \
               const _tname_ &, unsigned int, unsigned int,                    \
               const SubstructMatchParameters &, int) const) &                 \
               SubstructLibraryWrap::countMatches,                             \
           (python::arg("self"), python::arg("query"),                         \
            python::arg("startIdx"), python::arg("endIdx"),                    \
            python::arg("parameters"), python::arg("numThreads") = -1),        \
           "Get the matches for the query.\n\n"                                \
           " Arguments:\n"                                                     \
           "  - query:      substructure query\n"                              \
           "  - startIdx:   index to search from\n"                            \
           "  - endIdx:     index (non-inclusize) to search to\n"              \
           "  - numThreads: number of threads to use, -1 means all threads\n") \
      .def(                                                                    \
          "HasMatch",                                                          \
          (bool(SubstructLibraryWrap::*)(                                      \
              const _tname_ &, const SubstructMatchParameters &, int) const) & \
              SubstructLibraryWrap::hasMatch,                                  \
          (python::arg("self"), python::arg("query"),                          \
           python::arg("parameters"), python::arg("numThreads") = -1),         \
          "Get the matches for the query.\n\n"                                 \
          " Arguments:\n"                                                      \
          "  - query:      substructure query\n"                               \
          "  - numThreads: number of threads to use, -1 means all threads\n")  \
      .def("HasMatch",                                                         \
           (bool(SubstructLibraryWrap::*)(                                     \
               const _tname_ &, unsigned int, unsigned int,                    \
               const SubstructMatchParameters &, int) const) &                 \
               SubstructLibraryWrap::hasMatch,                                 \
           (python::arg("self"), python::arg("query"),                         \
            python::arg("startIdx"), python::arg("endIdx"),                    \
            python::arg("parameters"), python::arg("numThreads") = -1),        \
           "Get the matches for the query.\n\n"                                \
           " Arguments:\n"                                                     \
           "  - query:      substructure query\n"                              \
           "  - startIdx:   index to search from\n"                            \
           "  - endIdx:     index (non-inclusize) to search to\n"              \
           "  - numThreads: number of threads to use, -1 means all threads\n")

struct substructlibrary_wrapper {
  static void wrap() {
    python::class_<MolHolderBase, boost::shared_ptr<MolHolderBase>,
                   boost::noncopyable>("MolHolderBase", "", python::no_init)
        .def("__len__", &MolHolderBase::size, python::args("self"))
        .def("AddMol", &MolHolderBase::addMol, python::args("self", "m"),
             "Adds molecule to the molecule holder")
        .def("GetMol", &MolHolderBase::getMol, python::args("self", "arg1"),
             "Returns a particular molecule in the molecule holder\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which molecule to return\n\n"
             "    - sanitize: if sanitize is False, return the internal "
             "molecule state [default True]\n\n"
             "  NOTE: molecule indices start at 0\n")
        .def("__len__", &MolHolderBase::size, python::args("self"));

    python::class_<MolHolder, boost::shared_ptr<MolHolder>,
                   python::bases<MolHolderBase>>(
        "MolHolder", MolHolderDoc, python::init<>(python::args("self")));

    python::class_<CachedMolHolder, boost::shared_ptr<CachedMolHolder>,
                   python::bases<MolHolderBase>>(
        "CachedMolHolder", CachedMolHolderDoc,
        python::init<>(python::args("self")))
        .def("AddBinary", &CachedMolHolder::addBinary,
             (python::args("self", "pickle")),
             "Add a binary pickle to the molecule holder, no checking is done "
             "on the input data");

    python::class_<CachedSmilesMolHolder,
                   boost::shared_ptr<CachedSmilesMolHolder>,
                   python::bases<MolHolderBase>>(
        "CachedSmilesMolHolder", CachedSmilesMolHolderDoc,
        python::init<>(python::args("self")))
        .def("AddSmiles", &CachedSmilesMolHolder::addSmiles,
             (python::args("self", "smiles")),
             "Add a trusted smiles string to the molecule holder, no checking "
             "is done on the input data");

    python::class_<CachedTrustedSmilesMolHolder,
                   boost::shared_ptr<CachedTrustedSmilesMolHolder>,
                   python::bases<MolHolderBase>>(
        "CachedTrustedSmilesMolHolder", CachedTrustedSmilesMolHolderDoc,
        python::init<>(python::args("self")))
        .def("AddSmiles", &CachedTrustedSmilesMolHolder::addSmiles,
             (python::args("self", "smiles")),
             "Add a trusted smiles string to the molecule holder, no checking "
             "is done on the input data");

    python::class_<FPHolderBase, boost::shared_ptr<FPHolderBase>,
                   boost::noncopyable>("FPHolderBase", "", python::no_init)
        .def("__len__", &FPHolderBase::size, python::args("self"))

        .def("AddMol", &FPHolderBase::addMol, python::args("self", "m"),
             "Adds a molecule to the fingerprint database, returns the index "
             "of the new pattern")
        .def("AddFingerprint",
             (unsigned int (FPHolderBase::*)(const ExplicitBitVect &)) &
                 FPHolderBase::addFingerprint,
             python::args("self", "v"),
             "Adds a raw bit vector to the fingerprint database, returns the "
             "index of the supplied pattern")
        .def("GetFingerprint", &FPHolderBase::getFingerprint,
             python::return_value_policy<python::reference_existing_object>(),
             python::args("self", "idx"),
             "Return the bit vector at the specified index")
        .def("PassesFilter", &FPHolderBase::passesFilter,
             ((python::args("self"), python::args("idx")),
              python::args("query")),
             "Returns True if the specified index passes the filter supplied "
             "by the query bit vector")
        .def("MakeFingerprint", &FPHolderBase::makeFingerprint,
             ((python::arg("self"), python::arg("mol"))),
             python::return_value_policy<python::manage_new_object>(),
             "Compute the query bits for the holder");

    python::class_<PatternHolder, boost::shared_ptr<PatternHolder>,
                   python::bases<FPHolderBase>>(
        "PatternHolder", PatternHolderDoc, python::init<>(python::args("self")))
        .def(python::init<unsigned int>(python::args("self", "numBits")));

    python::class_<KeyHolderBase, boost::shared_ptr<KeyHolderBase>,
                   boost::noncopyable>("KeyHolderBase", "", python::no_init)
        .def("__len__", &KeyHolderBase::size, python::args("self"))

        .def("AddMol", &KeyHolderBase::addMol, python::args("self", "m"),
             "Adds a molecule to the fingerprint database, returns the index "
             "of the new pattern")
        .def("AddKey", &KeyHolderBase::addKey, python::args("self", "arg1"),
             "Add a key to the key holder, must be manually synced")
        .def("GetKey", &KeyHolderBase::getKey,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self", "arg1"),
             "Return the key at the specified index")
        .def("GetKeys", &KeyHolderBase::getKeys,
             python::args("self", "indices"),
             "Returns the keys for the given indices as return by GetMatches "
             "\n\n"
             "  ARGUMENTS:\n"
             "    - indices: The indices of the keys\n\n");

    python::class_<KeyFromPropHolder, boost::shared_ptr<KeyFromPropHolder>,
                   python::bases<KeyHolderBase>>(
        "KeyFromPropHolder", KeyHolderDoc, python::init<>(python::args("self")))
        .def(
            python::init<const std::string &>(python::args("self", "propname")))
        .def("GetPropName",
             (const std::string &(KeyFromPropHolder::*)() const) &
                 KeyFromPropHolder::getPropName,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"),
             "Return the key for the given molecule index");

    python::class_<TautomerPatternHolder,
                   boost::shared_ptr<TautomerPatternHolder>,
                   python::bases<FPHolderBase>>(
        "TautomerPatternHolder", TautomerPatternHolderDoc,
        python::init<>(python::args("self")))
        .def(python::init<unsigned int>(python::args("self", "numBits")));

    python::class_<SubstructLibraryWrap, boost::shared_ptr<SubstructLibraryWrap>>(
        "SubstructLibrary", SubstructLibraryDoc,
        python::init<>(python::args("self")))
        .def(python::init<boost::shared_ptr<MolHolderBase>>(
            python::args("self", "molecules")))
        .def(python::init<boost::shared_ptr<MolHolderBase>,
                          boost::shared_ptr<FPHolderBase>>(
            python::args("self", "molecules", "fingerprints")))
        .def(python::init<boost::shared_ptr<MolHolderBase>,
                          boost::shared_ptr<KeyHolderBase>>(
            python::args("self", "molecules", "keys")))
        .def(python::init<boost::shared_ptr<MolHolderBase>,
                          boost::shared_ptr<FPHolderBase>,
                          boost::shared_ptr<KeyHolderBase>>(
            python::args("self", "molecules", "fingerprints", "keys")))
        .def(python::init<std::string>(python::args("self", "pickle")))

        .def("GetMolHolder", &GetMolHolder, python::args("self"))
        .def("GetFpHolder", &GetFpHolder, python::args("self"))
        .def("GetKeyHolder", &GetKeyHolder, python::args("self"))

        .def("AddMol", &SubstructLibraryWrap::addMol,
             ((python::arg("self"), python::arg("mol"))),
             "Adds a molecule to the substruct library")

        // clang-format off
        LARGE_DEF(ROMol)
        LARGE_DEF(TautomerQuery)
        LARGE_DEF(MolBundle)
        LARGE_DEF(ExtendedQueryMol)
        // clang-format on

        .def("GetMol", &SubstructLibraryWrap::getMol,
             python::args("self", "idx"),
             "Returns a particular molecule in the molecule holder\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which molecule to return\n\n"
             "  NOTE: molecule indices start at 0\n")

        .def("SetSearchOrder", setSearchOrderHelper,
             python::args("self", "seq"),
             "Sets the search order for the library\n\n"
             "  ARGUMENTS:\n"
             "    - order: sequence of molecule indices\n\n"
             "  NOTE: molecule indices start at 0\n")
        .def("GetSearchOrder", getSearchOrderHelper, python::args("self"),
             "Returns the search order for the library\n\n"
             "  NOTE: molecule indices start at 0\n")

        .def("__len__", &SubstructLibraryWrap::size, python::args("self"))

        .def("ToStream", &toStream,
             (python::arg("self"), python::arg("stream")),
             "Serialize a substructure library to a python text stream.\n"
             "The stream can be a file in text mode or an io.StringIO type "
             "object\n\n"
             "  ARGUMENTS:\n"
             "    - stream: a text or text stream like object\n\n"
             "  >>> from rdkit.Chem import rdSubstructLibrary\n"
             "  >>> import io\n"
             "  >>> lib = rdSubstructLibrary.SubstructLibrary()\n"
             "  >>> stream = io.StringIO()\n"
             "  >>> lib.ToStream(stream)\n\n"
             "   or\n"
             "  >>> with open('rdkit.sslib', 'w') as stream:\n"
             "  ...  lib.ToStream(stream)\n")

        .def("InitFromStream", &initFromStream,
             (python::arg("self"), python::arg("stream")),
             "Deserialize a substructure library from a python bytes stream.\n"
             "Python doesn't allow seeking operations inside a unicode or "
             "string stream anymore\n"
             "so this requires opening a file in binary mode or using an "
             "io.ByteIO type object\n\n"
             "  ARGUMENTS:\n"
             "    - stream: a binary stream like object\n\n"
             "  SubstructLibrary.Serialize already writes a binary stream\n\n"
             "  >>> from rdkit.Chem import rdSubstructLibrary\n"
             "  >>> import io\n"
             "  >>> lib = rdSubstructLibrary.SubstructLibrary()\n"
             "  >>> stream = io.BytesIO( lib.Serialize() )\n"
             "  >>> lib.InitFromStream(stream)\n\n"
             "   remember to write to text and read from a binary stream\n"
             "  >>> with open('rdkit.sslib', 'w') as f: lib.ToStream(f)\n"
             "  >>> with open('rdkit.sslib', 'rb') as f: "
             "lib.InitFromStream(f)\n")

        .def("Serialize", &SubstructLibrary_Serialize, python::args("self"))
        // enable pickle support
        .def_pickle(substructlibrary_pickle_suite());

    python::def("SubstructLibraryCanSerialize", SubstructLibraryCanSerialize,
                "Returns True if the SubstructLibrary is serializable "
                "(requires boost serialization");

    python::def("AddPatterns",
                (void (*)(SubstructLibraryWrap &, int)) & addPatternsHelper,
                "Add pattern fingerprints to the given library, use "
                "numThreads=-1 to use all available cores",
                (python::arg("sslib"), python::arg("numThreads") = 1));

    python::def(
        "AddPatterns",
        (void (*)(SubstructLibraryWrap &, boost::shared_ptr<FPHolderBase>,
                  int)) &
            addPatternsHelper,
        "Add pattern fingerprints to the given library, use numThreads=-1 to "
        "use all available cores",
        (python::arg("sslib"), python::arg("patterns"),
         python::arg("numThreads") = 1));
  }
};
}  // namespace RDKit

void wrap_substructlibrary() { RDKit::substructlibrary_wrapper::wrap(); }
