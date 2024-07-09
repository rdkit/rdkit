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
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/FilterCatalog/FilterMatcherBase.h>
#include <GraphMol/FilterCatalog/FilterMatchers.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>
#include <GraphMol/RDKitBase.h>

namespace python = boost::python;

namespace RDKit {

struct filtercatalog_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const FilterCatalog &self) {
    std::string res;
    if (!FilterCatalogCanSerialize()) {
      throw_runtime_error("Pickling of FilterCatalog instances is not enabled");
    }
    res = self.Serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

template <typename T>
void python_to_vector(boost::python::object o, std::vector<T> &v) {
  python::stl_input_iterator<T> begin(o);
  python::stl_input_iterator<T> end;
  v.clear();
  v.insert(v.end(), begin, end);
}

void SetOffPatterns(ExclusionList &fc, boost::python::object list) {
  std::vector<FilterMatcherBase *> vect;
  // python_to_vector<FilterMatcherBase*>(list, vect);
  python::stl_input_iterator<FilterMatcherBase *> begin(list);
  python::stl_input_iterator<FilterMatcherBase *> end;

  std::vector<boost::shared_ptr<FilterMatcherBase>> temp;

  for (; begin != end; ++begin) {
    temp.push_back((*begin)->copy());
  }
  fc.setExclusionPatterns(temp);
}
/*
std::vector<const FilterCatalogEntry*> *GetMatches(FilterCatalog &fc, const
ROMol &mol) {
  std::vector<boost::shared_ptr<FilterCatalogEntry> > *result =        \
    new std::vector<const FilterCatalogEntry*>(fc.getMatches(mol));
  return result;
}
*/
std::vector<FilterMatch> FilterMatcherBaseGetMatches(FilterMatcherBase &fm,
                                                     const ROMol &mol) {
  std::vector<FilterMatch> matches;
  if (fm.getMatches(mol, matches)) {
    return matches;
  }
  return std::vector<FilterMatch>();
}

std::vector<FilterMatch> FilterCatalogEntryGetMatches(FilterCatalogEntry &fm,
                                                      const ROMol &mol) {
  std::vector<FilterMatch> matches;
  if (fm.getFilterMatches(mol, matches)) {
    return matches;
  }
  return std::vector<FilterMatch>();
}
/*
std::vector<MatchVectType> GetFilterMatchAtomPairs(FilterMatch &fm) {
  return fm.atomPairs;
}
*/
python::object FilterCatalogEntry_Serialize(const FilterCatalogEntry &cat) {
  std::string res = cat.Serialize();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object FilterCatalog_Serialize(const FilterCatalog &cat) {
  std::string res = cat.Serialize();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

int GetMatchVectItem(std::pair<int, int> &pair, size_t idx) {
  static const int def = 0xDEADBEEF;
  if (idx == 0) {
    return pair.first;
  } else if (idx == 1) {
    return pair.second;
  }
  PyErr_SetString(PyExc_IndexError, "Index out of bounds");
  python::throw_error_already_set();
  return def;
}

void filter_catalog_add_entry(FilterCatalog &catalog,
                              FilterCatalogEntry *entry) {
  // these are cheap to copy, so we do this to avoid memory management
  //  issues from python as the catalog will own the ptr
  catalog.addEntry(new FilterCatalogEntry(*entry));
}

class PythonFilterMatch : public FilterMatcherBase {
  PyObject *functor;
  bool incref;

 public:
  PythonFilterMatch(PyObject *self)
      : FilterMatcherBase("Python Filter Matcher"),
        functor(self),
        incref(false){};

  // ONLY CALLED FROM C++ from the copy operation
  PythonFilterMatch(const PythonFilterMatch &rhs)
      : FilterMatcherBase(rhs), functor(rhs.functor), incref(true) {
    python::incref(functor);
  }

  ~PythonFilterMatch() override {
    PyGILStateHolder h;
    if (incref) {
      python::decref(functor);
    }
  }
  bool isValid() const override {
    PyGILStateHolder h;
    return python::call_method<bool>(functor, "IsValid");
  }

  std::string getName() const override {
    PyGILStateHolder h;
    return python::call_method<std::string>(functor, "GetName");
  }

  bool getMatches(const ROMol &mol,
                  std::vector<FilterMatch> &matchVect) const override {
    PyGILStateHolder h;
    return python::call_method<bool>(functor, "GetMatches", boost::ref(mol),
                                     boost::ref(matchVect));
  }

  bool hasMatch(const ROMol &mol) const override {
    PyGILStateHolder h;
    return python::call_method<bool>(functor, "HasMatch", boost::ref(mol));
  }

  boost::shared_ptr<FilterMatcherBase> copy() const override {
    return boost::shared_ptr<FilterMatcherBase>(new PythonFilterMatch(*this));
  }
};

bool FilterCatalogRemoveEntry(FilterCatalog &fc,
                              const boost::python::object &obj) {
  if (PyLong_Check(obj.ptr())) {
    return fc.removeEntry(python::extract<unsigned int>(obj));
  }
  unsigned int idx =
      fc.getIdxForEntry(python::extract<FilterCatalogEntry *>(obj));
  return fc.removeEntry(idx);
}

const char *FilterMatchDoc =
    "Object that holds the result of running FilterMatcherBase::GetMatches\n\n"
    " - filterMatch holds the FilterMatchBase that triggered the match\n"
    " - atomParis holds the [ (query_atom_idx, target_atom_idx) ] pairs for "
    "the matches.\n"
    "\n\n"
    "Note that some matches may not have atom pairs (especially matches that "
    "use FilterMatchOps.Not";

const char *FilterMatcherBaseDoc =
    "Base class for matching molecules to filters.\n\n"
    " A FilterMatcherBase supplies the following API \n"
    " - IsValid() returns True if the matcher is valid for use, False "
    "otherwise\n"
    " - HasMatch(mol) returns True if the molecule matches the filter\n"
    " - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the "
    "FilterMatch data\n"
    "       that matches the molecule\n"
    "\n\nprint( FilterMatcherBase ) will print user-friendly information about "
    "the filter"
    "\nNote that a FilterMatcherBase can be combined from may "
    "FilterMatcherBases"
    "\nThis is why GetMatches can return multiple FilterMatcherBases.\n"
    ">>> from rdkit.Chem.FilterCatalog import *\n"
    ">>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', 0, 1)\n"
    ">>> oxygen_matcher = SmartsMatcher('Oxygen', '[#8]', 0, 1)\n"
    ">>> co_matcher = FilterMatchOps.Or(carbon_matcher, oxygen_matcher)\n"
    ">>> mol = Chem.MolFromSmiles('C')\n"
    ">>> matches = co_matcher.GetMatches(mol)\n"
    ">>> len(matches)\n"
    "1\n"
    ">>> print(matches[0].filterMatch)\n"
    "Carbon\n\n";

const char *SmartsMatcherDoc =
    "Smarts Matcher Filter\n"
    " basic constructors: \n"
    "   SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )\n"
    "   SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )\n\n"
    "  note: If the supplied smarts pattern is not valid, the IsValid() "
    "function will\n"
    "   return False\n"
    ">>> from rdkit.Chem.FilterCatalog import *\n"
    ">>> minCount, maxCount = 1,2\n"
    ">>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', minCount, maxCount)\n"
    ">>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CC')))\n"
    "True\n"
    ">>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))\n"
    "False\n"
    ">>> carbon_matcher.SetMinCount(2)\n"
    ">>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('C')))\n"
    "False\n"
    ">>> carbon_matcher.SetMaxCount(3)\n"
    ">>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))\n"
    "True\n"
    "\n";

const char *FilterHierarchyMatcherDoc =
    "Hierarchical Filter\n"
    " basic constructors: \n"
    "   FilterHierarchyMatcher( matcher )\n"
    "   where can be any FilterMatcherBase (SmartsMatcher, etc)\n"
    " FilterHierarchyMatcher's have children and can form matching\n"
    "  trees.  then GetFilterMatches is called, the most specific (\n"
    "  i.e. lowest node in a branch) is returned.\n\n"
    " n.b. A FilterHierarchicalMatcher of functional groups is returned\n"
    "  by calling GetFunctionalGroupHierarchy()\n\n"
    ">>> from rdkit.Chem import MolFromSmiles\n"
    ">>> from rdkit.Chem.FilterCatalog import *\n"
    ">>> functionalGroups = GetFunctionalGroupHierarchy()\n"
    ">>> [match.filterMatch.GetName() \n"
    "...     for match in functionalGroups.GetFilterMatches(\n"
    "...         MolFromSmiles('c1ccccc1Cl'))]\n"
    "['Halogen.Aromatic', 'Halogen.NotFluorine.Aromatic']\n"
    "\n";

const char *FilterCatalogEntryDoc =
    "FilterCatalogEntry\n"
    "A filter catalog entry is an entry in a filter catalog.\n"
    "Each filter is named and is used to flag a molecule usually for some\n"
    "undesirable property.\n\n"
    "For example, a PAINS (Pan Assay INterference) catalog entry be appear as\n"
    "follows:\n\n"
    ">>> from rdkit.Chem.FilterCatalog import *\n"
    ">>> params = FilterCatalogParams()\n"
    ">>> params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)\n"
    "True\n"
    ">>> catalog = FilterCatalog(params)\n"
    ">>> mol = "
    "Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/"
    "c1c(O)ccc2c1cccc2')\n"
    ">>> entry = catalog.GetFirstMatch(mol)\n"
    ">>> print (entry.GetProp('Scope'))\n"
    "PAINS filters (family A)\n"
    ">>> print (entry.GetDescription())\n"
    "hzone_phenol_A(479)\n"
    "\n\n";

python::dict GetFlattenedFunctionalGroupHierarchyHelper(bool normalize) {
  const std::map<std::string, ROMOL_SPTR> &flattened =
      GetFlattenedFunctionalGroupHierarchy(normalize);
  python::dict dict;
  for (const auto &it : flattened) {
    dict[it.first] = it.second;
  }
  return dict;
}

std::vector<std::vector<boost::shared_ptr<const FilterCatalogEntry>>>
RunFilterCatalogWrapper(const FilterCatalog &fc,
                        const std::vector<std::string> &smiles,
                        int numThreads) {
  NOGIL nogil;
  return RunFilterCatalog(fc, smiles, numThreads);
}

struct filtercat_wrapper {
  static void wrap() {
    python::class_<std::pair<int, int>>("IntPair")
        .def(python::init<const int &, const int &>(
            python::args("self", "query", "target")))
        .def_readwrite("query", &std::pair<int, int>::first)
        .def_readwrite("target", &std::pair<int, int>::second)
        .def("__getitem__", &GetMatchVectItem, python::args("self", "idx"));

    RegisterVectorConverter<std::pair<int, int>>("MatchTypeVect");

    python::class_<FilterMatch, boost::shared_ptr<FilterMatch>>(
        "FilterMatch", FilterMatchDoc,
        python::init<boost::shared_ptr<FilterMatcherBase>, MatchVectType>(
            python::args("self", "filter", "atomPairs")))
        .def_readonly("filterMatch", &FilterMatch::filterMatch)
        .def_readonly("atomPairs", &FilterMatch::atomPairs);

    RegisterVectorConverter<FilterMatch>("VectFilterMatch");

    python::class_<FilterMatcherBase,
                   boost::shared_ptr<FilterMatcherBase>, boost::noncopyable>(
        "FilterMatcherBase", FilterMatcherBaseDoc, python::no_init)
        .def("IsValid", &FilterMatcherBase::isValid, python::args("self"),
             "Return True if the filter matcher is valid, False otherwise")
        .def("HasMatch", &FilterMatcherBase::hasMatch,
             ((python::arg("self"), python::arg("mol"))),
             "Returns True if mol matches the filter")
        .def("GetMatches", &FilterMatcherBaseGetMatches,
             ((python::arg("self"), python::arg("mol"))),
             "Returns the list of matching subfilters mol matches any filter")

        .def("GetName", &FilterMatcherBase::getName, python::args("self"))
        .def("__str__", &FilterMatcherBase::getName, python::args("self"));

    python::register_ptr_to_python<boost::shared_ptr<FilterMatcherBase>>();

    python::class_<SmartsMatcher,
                   python::bases<FilterMatcherBase>>(
        "SmartsMatcher", SmartsMatcherDoc,
        python::init<const std::string &>(python::args("self", "name")))
        .def(python::init<const ROMol &>(python::args("self", "rhs"),
                                         "Construct from a molecule"))
        .def(python::init<const std::string &, const ROMol &, unsigned int,
                          unsigned int>(
            (python::arg("self"), python::arg("name"), python::arg("mol"),
             python::arg("minCount") = 1, python::arg("maxCount") = UINT_MAX),
            "Construct from a name, molecule, "
            "minimum and maximum count"))

        .def(python::init<const std::string &, const std::string &,
                          unsigned int, unsigned int>(
            (python::arg("self"), python::arg("name"), python::arg("smarts"),
             python::arg("minCount") = 1, python::arg("maxCount") = UINT_MAX),
            "Construct from a name,smarts pattern, minimum and "
            "maximum count"))

        .def("IsValid", &SmartsMatcher::isValid, python::args("self"),
             "Returns True if the SmartsMatcher is valid")
        .def(
            "SetPattern",
            (void(SmartsMatcher::*)(const ROMol &)) & SmartsMatcher::setPattern,
            python::args("self", "pat"),
            "Set the pattern molecule for the SmartsMatcher")
        .def("SetPattern",
             (void(SmartsMatcher::*)(const std::string &)) &
                 SmartsMatcher::setPattern,
             python::args("self", "pat"),
             "Set the smarts pattern for the Smarts Matcher (warning: "
             "MinimumCount is not reset)")
        .def("GetPattern", &SmartsMatcher::getPattern,
             python::return_value_policy<python::return_by_value>(),
             python::args("self"))
        .def("GetMinCount", &SmartsMatcher::getMinCount, python::args("self"),
             "Get the minimum times pattern must appear for the filter to "
             "match")
        .def("SetMinCount", &SmartsMatcher::setMinCount,
             ((python::arg("self"), python::arg("count"))),
             "Set the minimum times pattern must appear to match")

        .def("GetMaxCount", &SmartsMatcher::getMaxCount, python::args("self"),
             "Get the maximum times pattern can appear for the filter to match")
        .def(
            "SetMaxCount", &SmartsMatcher::setMaxCount,
            ((python::arg("self"), python::arg("count"))),
            "Set the maximum times pattern can appear for the filter to match");

    python::class_<ExclusionList,
                   python::bases<FilterMatcherBase>>(
        "ExclusionList", python::init<>(python::args("self")))
        .def("SetExclusionPatterns", &SetOffPatterns,
             python::args("self", "list"),
             "Set a list of FilterMatcherBases that should not appear in a "
             "molecule")
        .def("AddPattern", &ExclusionList::addPattern,
             python::args("self", "base"),
             "Add a FilterMatcherBase that should not appear in a molecule");

    python::class_<FilterHierarchyMatcher, FilterHierarchyMatcher *,
                   python::bases<FilterMatcherBase>>(
        "FilterHierarchyMatcher", FilterHierarchyMatcherDoc,
        python::init<>(python::args("self")))
        .def(python::init<const FilterMatcherBase &>(
            python::args("self", "matcher"), "Construct from a filtermatcher"))
        .def("SetPattern", &FilterHierarchyMatcher::setPattern,
             python::args("self", "matcher"),
             "Set the filtermatcher pattern for this node.  An empty node is "
             "considered "
             "a root node and passes along the matches to the children.")
        .def("AddChild", &FilterHierarchyMatcher::addChild,
             python::args("self", "hierarchy"),
             "Add a child node to this hierarchy.");

    python::register_ptr_to_python<boost::shared_ptr<FilterHierarchyMatcher>>();

    bool noproxy = true;
    RegisterVectorConverter<RDKit::ROMol *>("MolList", noproxy);

    python::class_<FilterCatalogEntry, FilterCatalogEntry *,
                   boost::shared_ptr<const FilterCatalogEntry>>(
        "FilterCatalogEntry", FilterCatalogEntryDoc,
        python::init<>(python::args("self")))
        .def(python::init<const std::string &, FilterMatcherBase &>(
            python::args("self", "name", "matcher")))
        .def("IsValid", &FilterCatalogEntry::isValid, python::args("self"))

        .def("GetDescription", &FilterCatalogEntry::getDescription,
             python::args("self"), "Get the description of the catalog entry")
        .def("SetDescription", &FilterCatalogEntry::setDescription,
             ((python::arg("self"), python::arg("description"))),
             "Set the description of the catalog entry")
        .def("GetFilterMatches", &FilterCatalogEntryGetMatches,
             (python::args("self", "mol")),
             "Retrieve the list of filters that match the molecule")
        .def("HasFilterMatch", &FilterCatalogEntry::hasFilterMatch,
             (python::args("self", "mol")),
             "Returns True if the catalog entry contains filters that match "
             "the molecule")

        .def("Serialize", &FilterCatalogEntry_Serialize, python::args("self"))
        .def("GetPropList", &FilterCatalogEntry::getPropList,
             python::args("self"))
        .def("SetProp",
             (void(FilterCatalogEntry::*)(const std::string &, std::string)) &
                 FilterCatalogEntry::setProp<std::string>,
             python::args("self", "key", "val"))
        .def("GetProp",
             (std::string(FilterCatalogEntry::*)(const std::string &) const) &
                 FilterCatalogEntry::getProp<std::string>,
             python::args("self", "key"))
        .def("ClearProp",
             (void(FilterCatalogEntry::*)(const std::string &)) &
                 FilterCatalogEntry::clearProp,
             python::args("self", "key"));

    python::register_ptr_to_python<boost::shared_ptr<FilterCatalogEntry>>();
    python::register_ptr_to_python<
        boost::shared_ptr<const FilterCatalogEntry>>();
    python::def(
        "GetFunctionalGroupHierarchy", GetFunctionalGroupHierarchy,
        "Returns the functional group hierarchy filter catalog",
        python::return_value_policy<python::reference_existing_object>());
    python::def(
        "GetFlattenedFunctionalGroupHierarchy",
        GetFlattenedFunctionalGroupHierarchyHelper,
        (python::args("normalized") = false),
        "Returns the flattened functional group hierarchy as a dictionary "
        " of name:ROMOL_SPTR substructure items");

    if (!is_python_converter_registered<
            boost::shared_ptr<const FilterCatalogEntry>>()) {
      python::register_ptr_to_python<
          boost::shared_ptr<const FilterCatalogEntry>>();
    }

    noproxy = true;
    RegisterVectorConverter<boost::shared_ptr<FilterCatalogEntry const>>(
        "FilterCatalogEntryList", noproxy);

    RegisterVectorConverter<
        std::vector<boost::shared_ptr<FilterCatalogEntry const>>>(
        "FilterCatalogListOfEntryList");

    {
      python::scope in_FilterCatalogParams =
          python::class_<FilterCatalogParams, FilterCatalogParams *>(
              "FilterCatalogParams", python::init<>(python::args("self")))
              .def(python::init<FilterCatalogParams::FilterCatalogs>(
                  python::args("self", "catalogs"),
                  "Construct from a FilterCatalogs identifier (i.e. "
                  "FilterCatalogParams.PAINS)"))
              .def("AddCatalog", &FilterCatalogParams::addCatalog,
                   python::args("self", "catalogs"));

      python::enum_<FilterCatalogParams::FilterCatalogs>("FilterCatalogs")
          .value("PAINS_A", FilterCatalogParams::PAINS_A)
          .value("PAINS_B", FilterCatalogParams::PAINS_B)
          .value("PAINS_C", FilterCatalogParams::PAINS_C)
          .value("PAINS", FilterCatalogParams::PAINS)
          .value("BRENK", FilterCatalogParams::BRENK)
          .value("NIH", FilterCatalogParams::NIH)
          .value("ZINC", FilterCatalogParams::ZINC)
          .value("CHEMBL_Glaxo", FilterCatalogParams::CHEMBL_Glaxo)
          .value("CHEMBL_Dundee", FilterCatalogParams::CHEMBL_Dundee)
          .value("CHEMBL_BMS", FilterCatalogParams::CHEMBL_BMS)
          .value("CHEMBL_SureChEMBL", FilterCatalogParams::CHEMBL_SureChEMBL)
          .value("CHEMBL_MLSMR", FilterCatalogParams::CHEMBL_MLSMR)
          .value("CHEMBL_Inpharmatica",
                 FilterCatalogParams::CHEMBL_Inpharmatica)
          .value("CHEMBL_LINT", FilterCatalogParams::CHEMBL_LINT)
          .value("CHEMBL", FilterCatalogParams::CHEMBL)
          .value("ALL", FilterCatalogParams::ALL);
    }

    python::class_<FilterCatalog>("FilterCatalog",
                                  python::init<>(python::args("self")))
        .def(python::init<const std::string &>(python::args("self", "binStr")))
        .def(python::init<const FilterCatalogParams &>(
            python::args("self", "params")))
        .def(python::init<FilterCatalogParams::FilterCatalogs>(
            python::args("self", "catalogs")))
        .def("Serialize", &FilterCatalog_Serialize, python::args("self"))
        .def("AddEntry", &filter_catalog_add_entry,
             (python::args("entry"), python::args("updateFPLength") = false),
             "Add a FilterCatalogEntry to the catalog")
        .def("RemoveEntry", &FilterCatalogRemoveEntry,
             python::args("self", "obj"),
             "Remove the given entry from the catalog")
        .def("GetNumEntries", &FilterCatalog::getNumEntries,
             python::args("self"),
             "Returns the number of entries in the catalog")
        .def("GetEntryWithIdx", &FilterCatalog::getEntry,
             ((python::arg("self"), python::arg("idx"))),
             "Return the FilterCatalogEntry at the specified index")
        .def("GetEntry", &FilterCatalog::getEntry,
             ((python::arg("self"), python::arg("idx"))),
             "Return the FilterCatalogEntry at the specified index")
        .def("HasMatch", &FilterCatalog::hasMatch,
             ((python::arg("self"), python::arg("mol"))),
             "Returns True if the catalog has an entry that matches mol")
        .def("GetFirstMatch", &FilterCatalog::getFirstMatch,
             ((python::arg("self"), python::arg("mol"))),
             "Return the first catalog entry that matches mol")
        .def("GetMatches", &FilterCatalog::getMatches,
             ((python::arg("self"), python::arg("mol"))),
             "Return all catalog entries that match mol")
        .def("GetFilterMatches", &FilterCatalog::getFilterMatches,
             ((python::arg("self"), python::arg("mol"))),
             "Return every matching filter from all catalog entries that match "
             "mol")
        // enable pickle support
        .def_pickle(filtercatalog_pickle_suite());

    python::class_<PythonFilterMatch, python::bases<FilterMatcherBase>>(
        "PythonFilterMatcher", python::init<PyObject *>(python::args("self")));

    python::def("FilterCatalogCanSerialize", FilterCatalogCanSerialize,
                "Returns True if the FilterCatalog is serializable "
                "(requires boost serialization");

    python::def("RunFilterCatalog", RunFilterCatalogWrapper,
                (python::arg("filterCatalog"), python::arg("smiles"),
                 python::arg("numThreads") = 1),
                "Run the filter catalog on the input list of smiles "
                "strings.\nUse numThreads=0 to use all available processors. "
                "Returns a vector of vectors.  For each input smiles, a vector "
                "of FilterCatalogEntry objects are "
                "returned for each matched filter.  If a molecule matches no "
                "filter, the vector will be empty. "
                "If a smiles string can't be parsed, a 'Bad smiles' entry is "
                "returned.");

    std::string nested_name = python::extract<std::string>(
        python::scope().attr("__name__") + ".FilterMatchOps");
    python::object nested_module(python::handle<>(
        python::borrowed(PyImport_AddModule(nested_name.c_str()))));
    python::scope().attr("FilterMatchOps") = nested_module;
    python::scope parent = nested_module;

    python::class_<FilterMatchOps::And,
                   python::bases<FilterMatcherBase>>(
        "And", python::init<FilterMatcherBase &, FilterMatcherBase &>(
                   python::args("self", "arg1", "arg2")));

    python::class_<FilterMatchOps::Or,
                   python::bases<FilterMatcherBase>>(
        "Or", python::init<FilterMatcherBase &, FilterMatcherBase &>(
                  python::args("self", "arg1", "arg2")));

    python::class_<FilterMatchOps::Not,
                   python::bases<FilterMatcherBase>>(
        "Not", python::init<FilterMatcherBase &>(python::args("self", "arg1")));
  };
};

}  // namespace RDKit

void wrap_filtercat() { RDKit::filtercat_wrapper::wrap(); }
