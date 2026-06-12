//
//  Copyright (C) 2015-2026 Novartis Institutes for BioMedical Research Inc.
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
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/map.h>
#include <nanobind/trampoline.h>
#include <RDBoost/boost_shared_ptr.h>

#include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/FilterCatalog/FilterMatcherBase.h>
#include <GraphMol/FilterCatalog/FilterMatchers.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
std::string serializeFilterCatalog(const RDKit::FilterCatalog &catalog) {
  return catalog.Serialize();
}
}  // namespace

namespace RDKit {

// Trampoline class to allow Python subclassing of FilterMatcherBase.
// We hold a strong nb::object reference to the Python side to keep it alive,
// which is necessary when the C++ copy() method creates a copied instance
// (e.g. FilterCatalogEntry stores a copy).
struct FilterMatcherBaseTrampoline : FilterMatcherBase {
  NB_TRAMPOLINE(FilterMatcherBase, 4);
  // Strong reference to the Python object — keeps it alive when the C++
  // side is held in a boost::shared_ptr without a Python wrapper.
  mutable nb::object d_pyObject;

  bool isValid() const override {
    nb::gil_scoped_acquire gil;
    nb::object pyObj = d_pyObject.is_valid()
                           ? d_pyObject
                           : nb::borrow<nb::object>(nb_trampoline.base());
    return nb::cast<bool>(pyObj.attr("IsValid")());
  }

  std::string getName() const override {
    nb::gil_scoped_acquire gil;
    nb::object pyObj = d_pyObject.is_valid()
                           ? d_pyObject
                           : nb::borrow<nb::object>(nb_trampoline.base());
    return nb::cast<std::string>(pyObj.attr("GetName")());
  }

  bool getMatches(const ROMol &mol,
                  std::vector<FilterMatch> &matchVect) const override {
    // Call Python's GetMatches(mol, vect) where vect is a Python list.
    // The Python implementation appends FilterMatch objects to vect.
    // We then copy those back into the C++ matchVect.
    nb::gil_scoped_acquire gil;
    nb::object pyObj = d_pyObject.is_valid()
                           ? d_pyObject
                           : nb::borrow<nb::object>(nb_trampoline.base());
    nb::list pyVect;
    nb::object result = pyObj.attr("GetMatches")(mol, pyVect);
    bool matched = nb::cast<bool>(result);
    if (matched) {
      for (auto item : pyVect) {
        matchVect.push_back(nb::cast<FilterMatch>(item));
      }
    }
    return matched;
  }

  bool hasMatch(const ROMol &mol) const override {
    nb::gil_scoped_acquire gil;
    nb::object pyObj = d_pyObject.is_valid()
                           ? d_pyObject
                           : nb::borrow<nb::object>(nb_trampoline.base());
    return nb::cast<bool>(pyObj.attr("HasMatch")(mol));
  }

  boost::shared_ptr<FilterMatcherBase> copy() const override {
    nb::gil_scoped_acquire gil;
    // Capture the Python object so the copy keeps it alive.
    nb::object pyObj = d_pyObject.is_valid()
                           ? d_pyObject
                           : nb::borrow<nb::object>(nb_trampoline.base());
    auto *copy = new FilterMatcherBaseTrampoline(*this);
    copy->d_pyObject = pyObj;
    return boost::shared_ptr<FilterMatcherBase>(copy);
  }
};

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

void SetOffPatterns(ExclusionList &fc, nb::object list) {
  std::vector<boost::shared_ptr<FilterMatcherBase>> temp;
  for (auto item : list) {
    FilterMatcherBase *matcher = nb::cast<FilterMatcherBase *>(item);
    temp.push_back(matcher->copy());
  }
  fc.setExclusionPatterns(temp);
}

void filter_catalog_add_entry(FilterCatalog &catalog,
                              FilterCatalogEntry *entry) {
  catalog.addEntry(new FilterCatalogEntry(*entry));
}

bool FilterCatalogRemoveEntry(FilterCatalog &fc, nb::object obj) {
  if (PyLong_Check(obj.ptr())) {
    return fc.removeEntry(nb::cast<unsigned int>(obj));
  }
  unsigned int idx = fc.getIdxForEntry(nb::cast<FilterCatalogEntry *>(obj));
  return fc.removeEntry(idx);
}

nb::dict GetFlattenedFunctionalGroupHierarchyHelper(bool normalize) {
  const std::map<std::string, ROMOL_SPTR> &flattened =
      GetFlattenedFunctionalGroupHierarchy(normalize);
  nb::dict dict;
  for (const auto &it : flattened) {
    dict[it.first.c_str()] = it.second;
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

void wrap_filtercat(nb::module_ &m) {
  // MatchTypeVect is a Python list (std::vector<std::pair<int,int>>
  // auto-converts).
  m.attr("MatchTypeVect") = nb::module_::import_("builtins").attr("list");
  // IntPair is a factory that takes two ints and returns a (int, int) tuple.
  m.def(
      "IntPair", [](int a, int b) { return std::make_pair(a, b); }, "a"_a,
      "b"_a, "Create an integer pair (tuple) for use in MatchTypeVect");

  nb::class_<FilterMatch>(
      m, "FilterMatch",
      R"DOC(Object that holds the result of running FilterMatcherBase::GetMatches

 - filterMatch holds the FilterMatchBase that triggered the match
 - atomPairs holds the [ (query_atom_idx, target_atom_idx) ] pairs for the matches.

Note that some matches may not have atom pairs (especially matches that
use FilterMatchOps.Not)DOC")
      .def(nb::init<boost::shared_ptr<FilterMatcherBase>, MatchVectType>(),
           "filter"_a, "atomPairs"_a)
      .def_ro("filterMatch", &FilterMatch::filterMatch)
      .def_ro("atomPairs", &FilterMatch::atomPairs);

  nb::class_<FilterMatcherBase, FilterMatcherBaseTrampoline>(
      m, "FilterMatcher",
      R"DOC(Base class for matching molecules to filters.

 A FilterMatcherBase supplies the following API
 - IsValid() returns True if the matcher is valid for use, False otherwise
 - HasMatch(mol) returns True if the molecule matches the filter
 - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the FilterMatch data
       that matches the molecule

print( FilterMatcherBase ) will print user-friendly information about the filter
Note that a FilterMatcherBase can be combined from many FilterMatcherBases
This is why GetMatches can return multiple FilterMatcherBases.
>>> from rdkit.Chem.FilterCatalog import *
>>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', 0, 1)
>>> oxygen_matcher = SmartsMatcher('Oxygen', '[#8]', 0, 1)
>>> co_matcher = FilterMatchOps.Or(carbon_matcher, oxygen_matcher)
>>> mol = Chem.MolFromSmiles('C')
>>> matches = co_matcher.GetMatches(mol)
>>> len(matches)
1
>>> print(matches[0].filterMatch)
Carbon
)DOC")
      .def(nb::init<const std::string &>(), "name"_a)
      .def("IsValid", &FilterMatcherBase::isValid,
           "Return True if the filter matcher is valid, False otherwise")
      .def("HasMatch", &FilterMatcherBase::hasMatch, "mol"_a,
           "Returns True if mol matches the filter")
      .def("GetMatches", &FilterMatcherBaseGetMatches, "mol"_a,
           "Returns the list of matching subfilters mol matches any filter")
      .def("GetName", &FilterMatcherBase::getName)
      .def("__str__", &FilterMatcherBase::getName);

  // Also expose FilterMatcherBase and PythonFilterMatcher under their original
  // names for compatibility
  m.attr("FilterMatcherBase") = m.attr("FilterMatcher");
  m.attr("PythonFilterMatcher") = m.attr("FilterMatcher");

  nb::class_<SmartsMatcher, FilterMatcherBase>(m, "SmartsMatcher",
                                               R"DOC(Smarts Matcher Filter
 basic constructors:
   SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )
   SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )

  note: If the supplied smarts pattern is not valid, the IsValid() function will
   return False
>>> from rdkit.Chem.FilterCatalog import *
>>> minCount, maxCount = 1,2
>>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', minCount, maxCount)
>>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CC')))
True
>>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
False
>>> carbon_matcher.SetMinCount(2)
>>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('C')))
False
>>> carbon_matcher.SetMaxCount(3)
>>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
True
)DOC")
      .def(nb::init<const std::string &>(), "name"_a)
      .def(nb::init<const ROMol &>(), "rhs"_a, "Construct from a molecule")
      .def(nb::init<const std::string &, const ROMol &, unsigned int,
                    unsigned int>(),
           "name"_a, "mol"_a, "minCount"_a = 1, "maxCount"_a = UINT_MAX,
           "Construct from a name, molecule, minimum and maximum count")
      .def(nb::init<const std::string &, const std::string &, unsigned int,
                    unsigned int>(),
           "name"_a, "smarts"_a, "minCount"_a = 1, "maxCount"_a = UINT_MAX,
           "Construct from a name, smarts pattern, minimum and maximum count")
      .def("IsValid", &SmartsMatcher::isValid,
           "Returns True if the SmartsMatcher is valid")
      .def("SetPattern",
           nb::overload_cast<const ROMol &>(&SmartsMatcher::setPattern),
           "pat"_a, "Set the pattern molecule for the SmartsMatcher")
      .def("SetPattern",
           nb::overload_cast<const std::string &>(&SmartsMatcher::setPattern),
           "pat"_a,
           "Set the smarts pattern for the Smarts Matcher (warning: "
           "MinimumCount is not reset)")
      .def("GetPattern", &SmartsMatcher::getPattern)
      .def("GetMinCount", &SmartsMatcher::getMinCount,
           "Get the minimum times pattern must appear for the filter to match")
      .def("SetMinCount", &SmartsMatcher::setMinCount, "count"_a,
           "Set the minimum times pattern must appear to match")
      .def("GetMaxCount", &SmartsMatcher::getMaxCount,
           "Get the maximum times pattern can appear for the filter to match")
      .def("SetMaxCount", &SmartsMatcher::setMaxCount, "count"_a,
           "Set the maximum times pattern can appear for the filter to match");

  nb::class_<ExclusionList, FilterMatcherBase>(m, "ExclusionList")
      .def(nb::init<>())
      .def("SetExclusionPatterns", &SetOffPatterns, "list"_a,
           "Set a list of FilterMatcherBases that should not appear in a "
           "molecule")
      .def("AddPattern", &ExclusionList::addPattern, "base"_a,
           "Add a FilterMatcherBase that should not appear in a molecule");

  nb::class_<FilterHierarchyMatcher, FilterMatcherBase>(
      m, "FilterHierarchyMatcher",
      R"DOC(Hierarchical Filter
 basic constructors:
   FilterHierarchyMatcher( matcher )
   where can be any FilterMatcherBase (SmartsMatcher, etc)
 FilterHierarchyMatcher's have children and can form matching
  trees.  When GetFilterMatches is called, the most specific (
  i.e. lowest node in a branch) is returned.

 n.b. A FilterHierarchicalMatcher of functional groups is returned
  by calling GetFunctionalGroupHierarchy()

>>> from rdkit.Chem import MolFromSmiles
>>> from rdkit.Chem.FilterCatalog import *
>>> functionalGroups = GetFunctionalGroupHierarchy()
>>> [match.filterMatch.GetName()
...     for match in functionalGroups.GetFilterMatches(
...         MolFromSmiles('c1ccccc1Cl'))]
['Halogen.Aromatic', 'Halogen.NotFluorine.Aromatic']
)DOC")
      .def(nb::init<>())
      .def(nb::init<const FilterMatcherBase &>(), "matcher"_a,
           "Construct from a filtermatcher")
      .def("SetPattern", &FilterHierarchyMatcher::setPattern, "matcher"_a,
           R"DOC(Set the filtermatcher pattern for this node. An empty node is
considered a root node and passes along the matches to the children.)DOC")
      .def("AddChild", &FilterHierarchyMatcher::addChild, "hierarchy"_a,
           "Add a child node to this hierarchy.");

  nb::class_<FilterCatalogEntry>(m, "FilterCatalogEntry",
                                 R"DOC(FilterCatalogEntry
A filter catalog entry is an entry in a filter catalog.
Each filter is named and is used to flag a molecule usually for some
undesirable property.

For example, a PAINS (Pan Assay INterference) catalog entry be appear as
follows:

>>> from rdkit.Chem.FilterCatalog import *
>>> params = FilterCatalogParams()
>>> params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
True
>>> catalog = FilterCatalog(params)
>>> mol = Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2')
>>> entry = catalog.GetFirstMatch(mol)
>>> print (entry.GetProp('Scope'))
PAINS filters (family A)
>>> print (entry.GetDescription())
hzone_phenol_A(479)
)DOC")
      .def(nb::init<>())
      .def(nb::init<const std::string &, FilterMatcherBase &>(), "name"_a,
           "matcher"_a)
      .def("IsValid", &FilterCatalogEntry::isValid)
      .def("GetDescription", &FilterCatalogEntry::getDescription,
           "Get the description of the catalog entry")
      .def("SetDescription", &FilterCatalogEntry::setDescription,
           "description"_a, "Set the description of the catalog entry")
      .def("GetFilterMatches", &FilterCatalogEntryGetMatches, "mol"_a,
           "Retrieve the list of filters that match the molecule")
      .def("HasFilterMatch", &FilterCatalogEntry::hasFilterMatch, "mol"_a,
           "Returns True if the catalog entry contains filters that match "
           "the molecule")
      .def("Serialize",
           [](const FilterCatalogEntry &cat) {
             const auto res = cat.Serialize();
             return nb::bytes(res.c_str(), res.size());
           })
      .def("GetPropList", &FilterCatalogEntry::getPropList)
      .def(
          "SetProp",
          [](FilterCatalogEntry &entry, const std::string &key,
             const std::string &val) { entry.setProp<std::string>(key, val); },
          "key"_a, "val"_a)
      .def(
          "GetProp",
          [](const FilterCatalogEntry &entry, const std::string &key) {
            return entry.getProp<std::string>(key);
          },
          "key"_a)
      .def(
          "ClearProp",
          [](FilterCatalogEntry &entry, const std::string &key) {
            entry.clearProp(key);
          },
          "key"_a);

  m.def("GetFunctionalGroupHierarchy", &GetFunctionalGroupHierarchy,
        "Returns the functional group hierarchy filter catalog",
        nb::rv_policy::reference);
  m.def("GetFlattenedFunctionalGroupHierarchy",
        &GetFlattenedFunctionalGroupHierarchyHelper, "normalized"_a = false,
        "Returns the flattened functional group hierarchy as a dictionary "
        " of name:ROMOL_SPTR substructure items");

  {
    nb::class_<FilterCatalogParams> params_class(m, "FilterCatalogParams");
    params_class.def(nb::init<>())
        .def(nb::init<FilterCatalogParams::FilterCatalogs>(), "catalogs"_a,
             "Construct from a FilterCatalogs identifier (i.e. "
             "FilterCatalogParams.PAINS)")
        .def("AddCatalog", &FilterCatalogParams::addCatalog, "catalogs"_a);

    nb::enum_<FilterCatalogParams::FilterCatalogs>(
        params_class, "FilterCatalogs", nb::is_arithmetic())
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
        .value("CHEMBL_Inpharmatica", FilterCatalogParams::CHEMBL_Inpharmatica)
        .value("CHEMBL_LINT", FilterCatalogParams::CHEMBL_LINT)
        .value("CHEMBL", FilterCatalogParams::CHEMBL)
        .value("ALL", FilterCatalogParams::ALL);
  }

  nb::class_<FilterCatalog>(m, "FilterCatalog")
      .def(nb::new_([]() { return new FilterCatalog(); }))
      .def(nb::new_([](nb::bytes pkl) {
             return new FilterCatalog(std::string(
                 static_cast<const char *>(pkl.data()), pkl.size()));
           }),
           "pickle"_a)
      .def(nb::new_([](const FilterCatalogParams &ps) {
             return new FilterCatalog(ps);
           }),
           "params"_a)
      .def(nb::new_([](FilterCatalogParams::FilterCatalogs fcs) {
             return new FilterCatalog(fcs);
           }),
           "catalogs"_a)
      .def("Serialize",
           [](const FilterCatalog &cat) {
             if (!FilterCatalogCanSerialize()) {
               throw std::runtime_error(
                   "Pickling of FilterCatalog instances is not enabled");
             }
             const auto res = cat.Serialize();
             return nb::bytes(res.c_str(), res.size());
           })
      .def("AddEntry", &filter_catalog_add_entry, "entry"_a,
           "Add a FilterCatalogEntry to the catalog")
      .def("RemoveEntry", &FilterCatalogRemoveEntry, "obj"_a,
           "Remove the given entry from the catalog")
      .def("GetNumEntries", &FilterCatalog::getNumEntries,
           "Returns the number of entries in the catalog")
      .def("GetEntryWithIdx", &FilterCatalog::getEntry, "idx"_a,
           "Return the FilterCatalogEntry at the specified index")
      .def("GetEntry", &FilterCatalog::getEntry, "idx"_a,
           "Return the FilterCatalogEntry at the specified index")
      .def("HasMatch", &FilterCatalog::hasMatch, "mol"_a,
           "Returns True if the catalog has an entry that matches mol")
      .def("GetFirstMatch", &FilterCatalog::getFirstMatch, "mol"_a,
           "Return the first catalog entry that matches mol")
      .def("GetMatches", &FilterCatalog::getMatches, "mol"_a,
           "Return all catalog entries that match mol")
      .def("GetFilterMatches", &FilterCatalog::getFilterMatches, "mol"_a,
           "Return every matching filter from all catalog entries that match "
           "mol")
      .def("__setstate__", setObjectState<FilterCatalog>)
      .def("__getstate__",
           getObjectState<FilterCatalog, serializeFilterCatalog>);

  m.def("FilterCatalogCanSerialize", FilterCatalogCanSerialize,
        "Returns True if the FilterCatalog is serializable "
        "(requires boost serialization)");

  m.def("RunFilterCatalog", &RunFilterCatalogWrapper, "filterCatalog"_a,
        "smiles"_a, "numThreads"_a = 1,
        R"DOC(Run the filter catalog on the input list of smiles strings.
Use numThreads=0 to use all available processors.
Returns a vector of vectors. For each input smiles, a vector of
FilterCatalogEntry objects are returned for each matched filter. If a
molecule matches no filter, the vector will be empty. If a smiles string
can't be parsed, a 'Bad smiles' entry is returned.)DOC");

  // Create the FilterMatchOps submodule
  nb::module_ ops_module(nb::steal<nb::module_>(
      PyImport_AddModule("rdkit.Chem.rdfiltercatalog.FilterMatchOps")));
  m.attr("FilterMatchOps") = ops_module;

  nb::class_<FilterMatchOps::And, FilterMatcherBase>(ops_module, "And")
      .def(nb::init<FilterMatcherBase &, FilterMatcherBase &>(), "arg1"_a,
           "arg2"_a);

  nb::class_<FilterMatchOps::Or, FilterMatcherBase>(ops_module, "Or")
      .def(nb::init<FilterMatcherBase &, FilterMatcherBase &>(), "arg1"_a,
           "arg2"_a);

  nb::class_<FilterMatchOps::Not, FilterMatcherBase>(ops_module, "Not")
      .def(nb::init<FilterMatcherBase &>(), "arg1"_a);
}

}  // namespace RDKit
