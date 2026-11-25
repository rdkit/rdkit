//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/optional.h>
#include <nanobind/trampoline.h>
#include <nanobind/ndarray.h>

#include <RDBoost/Wrap_nb.h>
#include "../NormalizedDescriptors.h"

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit::NormalizedDescriptors;

namespace {
CumulativeHistogram *makeCumulativeHistogram(const nb::object &edges,
                                             const nb::object &fractions) {
  std::vector<double> edgev, fracv;
  pythonObjectToVect<double>(edges, edgev);
  pythonObjectToVect<double>(fractions, fracv);
  return new CumulativeHistogram(std::move(edgev), std::move(fracv));
}

nb::tuple toTuple(const std::vector<double> &vals) {
  nb::list res;
  for (auto v : vals) {
    res.append(v);
  }
  return nb::tuple(res);
}

nb::tuple getEdges(const CumulativeHistogram &table) {
  return toTuple(table.getEdges());
}

nb::tuple getFractions(const CumulativeHistogram &table) {
  return toTuple(table.getFractions());
}

nb::list getNames(const HistogramSet &tables) {
  nb::list res;
  for (const auto &name : tables.getNames()) {
    res.append(name);
  }
  return res;
}

void loadFromString(HistogramSet &tables, const std::string &text) {
  std::istringstream inStream(text);
  tables.loadFromStream(inStream);
}

std::string toString(const HistogramSet &tables) {
  std::ostringstream outStream;
  tables.writeToStream(outStream);
  return outStream.str();
}

nb::list toList(const std::vector<double> &vals) {
  nb::list res;
  for (auto v : vals) {
    res.append(v);
  }
  return res;
}

nb::list calcNormalizedHelper(const RDKit::ROMol &mol, nb::object tables) {
  if (tables.is_none()) {
    return toList(calcNormalizedDescriptors(mol));
  }
  return toList(
      calcNormalizedDescriptors(mol, nb::cast<const HistogramSet &>(tables)));
}

nb::list calcDescriptorValuesHelper(const RDKit::ROMol &mol) {
  return toList(calcDescriptorValues(mol));
}

nb::tuple getNormalizedDescriptorNamesHelper() {
  nb::list res;
  for (const auto &name : getNormalizedDescriptorNames()) {
    res.append(name);
  }
  return nb::tuple(res);
}

NormalizedProperties *makeNormalizedProperties(nb::object names,
                                               nb::object tables) {
  std::shared_ptr<const HistogramSet> tableSet;
  if (tables.is_none()) {
    tableSet.reset(&getDefaultTables(), [](const HistogramSet *) {});
  } else {
    tableSet = std::make_shared<const HistogramSet>(
        nb::cast<const HistogramSet &>(tables));
  }
  if (names.is_none()) {
    return new NormalizedProperties(getNormalizedDescriptorNames(), tableSet);
  }
  std::vector<std::string> nameVect;
  pythonObjectToVect<std::string>(names, nameVect);
  return new NormalizedProperties(nameVect, tableSet);
}

nb::list computeRawHelper(const NormalizedProperties &props,
                          const RDKit::ROMol &mol) {
  return toList(props.computeRawProperties(mol));
}

double normalizeWithDefaults(const std::string &name, double value) {
  return getDefaultTables().normalize(name, value);
}
}  // namespace

NB_MODULE(rdNormalizedDescriptors, m) {
  m.doc() =
      "Module containing histogram normalization of molecular descriptors, a "
      "port of the RDKit2DHistogramNormalized descriptors from "
      "descriptastorus";

  nb::class_<CumulativeHistogram>(
      m, "CumulativeHistogram",
      "A cumulative histogram for one descriptor: the left edge of each bin "
      "and the fraction of the reference set in that bin or below it.\n"
      "A value normalizes to the fraction of the first bin whose edge is >= "
      "the value, or 1.0 past the last edge. NaN normalizes to 0.0.")
      .def(nb::init<const std::vector<double>, std::vector<double>>(),
           "edges"_a, "fractions"_a)
      .def("Normalize", &CumulativeHistogram::normalize, "value"_a,
           "returns the normalized value")
      .def("__call__", &CumulativeHistogram::normalize, "value"_a)
      .def("GetEdges", getEdges)
      .def("GetFractions", getFractions);

  nb::class_<HistogramSet>(
      m, "HistogramSet",
      "A collection of histogram tables keyed by descriptor name")
      .def(nb::init<>())
      .def("AddTable", &HistogramSet::addTable, "name"_a, "table"_a,
           "adds (or replaces) the table for a descriptor")
      .def("HasTable", &HistogramSet::hasTable, ("name"_a))
      .def("GetTable",
           (const CumulativeHistogram &(HistogramSet::*)(size_t) const) &
               HistogramSet::getTable,
           "idx"_a,
           "returns the table with an index, raises IndexError if out of "
           "range")
      .def("GetTable",
           (const CumulativeHistogram &(HistogramSet::*)(const std::string &)
                const) &
               HistogramSet::getTable,
           "name"_a,
           "returns the table for a descriptor, raises KeyError if missing")
      .def("GetNames", getNames, "returns the descriptor names, in index order")
      .def("GetTableIndex", &HistogramSet::getTableIndex, "name"_a,
           "returns the index of the table for a descriptor, -1 if missing")
      .def("Normalize",
           (double (HistogramSet::*)(size_t, double) const) &
               HistogramSet::normalize,
           "idx"_a, "value"_a,
           "normalizes a value with the table with an index")
      .def("Normalize",
           (double (HistogramSet::*)(const std::string &, double) const) &
               HistogramSet::normalize,
           "name"_a, "value"_a,
           "normalizes a descriptor value, returns 0.0 for descriptors "
           "without a table")
      .def("LoadFromFile", &HistogramSet::loadFromFile, "fileName"_a,
           "adds the tables in a file to this set")
      .def("LoadFromString", loadFromString, "text"_a,
           "adds the tables in a string to this set")
      .def("ToString", toString,
           "returns the tables in the text format read by LoadFromString")
      .def("__len__", &HistogramSet::size);

   nb::module_::import_("rdkit.Chem.rdMolDescriptors");
   
  nb::class_<NormalizedProperties, RDKit::Descriptors::Properties>(
      m, "NormalizedProperties",
      "Computes descriptors with the rdMolDescriptors.Properties registry and "
      "normalizes them with histogram tables.\n"
      "names defaults to GetNormalizedDescriptorNames() and tables to "
      "GetDefaultTables(). ComputeProperties returns normalized values;")
    .def(nb::init<>())
    .def("__init__", [](NormalizedProperties *self) {
      new (self) NormalizedProperties(getNormalizedDescriptorNames(), getDefaultTables());
    })
    .def("__init__", [](NormalizedProperties *self,
			const std::vector<std::string> &names) {
      new (self) NormalizedProperties(names, getDefaultTables());
    }, "names"_a)
    .def("__init__", [](NormalizedProperties *self,
			const std::vector<std::string> &names, const HistogramSet &tables) {
      new (self) NormalizedProperties(names, tables);
    },
	 "names"_a, "tables"_a = nb::none())
    .def("ComputeRawProperties", computeRawHelper, "mol"_a,
	 "returns the unnormalized descriptor values; failures are nan")
    .def("GetTables", &NormalizedProperties::getTables);
  
  m.def("GetPropertyName", getPropertyName, "name"_a,
              "returns the name of the rdMolDescriptors property that computes "
              "a descriptor, e.g. 'exactmw' for 'ExactMolWt'");

  m.def("GetDefaultTablePath", getDefaultTablePath,
        "returns the path of the default histogram tables");
  m.def("GetDefaultTables", getDefaultTables,
        "returns the default histogram tables (descriptastorus' "
        "RDKit2DHistogramNormalized reference distributions)");
  m.def("GetNormalizedDescriptorNames", getNormalizedDescriptorNamesHelper,
        "returns the names of the descriptors calculated by "
        "CalcNormalizedDescriptors, in the order they are returned "
        "(the order of rdkit.Chem.Descriptors._descList)");
  m.def("CalcNormalizedDescriptors", calcNormalizedHelper, "mol"_a,
        "tables"_a = nb::none(),
        "calculates the descriptors named by "
        "GetNormalizedDescriptorNames and normalizes them "
        "with tables (the default tables when None). Descriptors that "
        "cannot be calculated or have no table are 0.0");
  m.def("CalcDescriptorValues", calcDescriptorValuesHelper, "mol"_a,
        "calculates the raw (unnormalized) values of the descriptors "
        "named by GetNormalizedDescriptorNames; failures are nan");
  m.def("NormalizeDescriptor", normalizeWithDefaults, "name"_a, "value"_a,
        "normalizes a descriptor value with the default tables, returns "
        "0.0 for descriptors without a table");
}
