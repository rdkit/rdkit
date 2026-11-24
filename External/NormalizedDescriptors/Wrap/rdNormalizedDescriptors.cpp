//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/Property.h>

#include <memory>
#include "../NormalizedDescriptors.h"

namespace python = boost::python;
using namespace RDKit::NormalizedDescriptors;

namespace {
CumulativeHistogram *makeCumulativeHistogram(const python::object &edges,
                                             const python::object &fractions) {
  std::vector<double> edgev, fracv;
  pythonObjectToVect<double>(edges, edgev);
  pythonObjectToVect<double>(fractions, fracv);
  return new CumulativeHistogram(std::move(edgev), std::move(fracv));
}

python::tuple toTuple(const std::vector<double> &vals) {
  python::list res;
  for (auto v : vals) {
    res.append(v);
  }
  return python::tuple(res);
}

python::tuple getEdges(const CumulativeHistogram &table) {
  return toTuple(table.getEdges());
}

python::tuple getFractions(const CumulativeHistogram &table) {
  return toTuple(table.getFractions());
}

python::list getNames(const HistogramSet &tables) {
  python::list res;
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

python::list toList(const std::vector<double> &vals) {
  python::list res;
  for (auto v : vals) {
    res.append(v);
  }
  return res;
}

python::list calcNormalizedHelper(const RDKit::ROMol &mol,
                                  python::object tables) {
  if (tables.is_none()) {
    return toList(calcNormalizedDescriptors(mol));
  }
  python::extract<const HistogramSet &> tableSet(tables);
  return toList(calcNormalizedDescriptors(mol, tableSet()));
}

python::list calcDescriptorValuesHelper(const RDKit::ROMol &mol) {
  return toList(calcDescriptorValues(mol));
}

python::tuple getNormalizedDescriptorNamesHelper() {
  python::list res;
  for (const auto &name : getNormalizedDescriptorNames()) {
    res.append(name);
  }
  return python::tuple(res);
}

NormalizedProperties *makeNormalizedProperties(python::object names,
                                               python::object tables) {
  std::shared_ptr<const HistogramSet> tableSet;
  if (tables.is_none()) {
    tableSet.reset(&getDefaultTables(), [](const HistogramSet *) {});
  } else {
    tableSet = std::make_shared<const HistogramSet>(
        python::extract<const HistogramSet &>(tables)());
  }
  if (names.is_none()) {
    return new NormalizedProperties(getNormalizedDescriptorNames(), tableSet);
  }
  std::vector<std::string> nameVect;
  pythonObjectToVect<std::string>(names, nameVect);
  return new NormalizedProperties(nameVect, tableSet);
}

python::tuple getDescriptorNames(const NormalizedProperties &props) {
  python::list res;
  for (const auto &name : props.getDescriptorNames()) {
    res.append(name);
  }
  return python::tuple(res);
}

python::list computeRawHelper(const NormalizedProperties &props,
                              const RDKit::ROMol &mol) {
  return toList(props.computeRawProperties(mol));
}

double normalizeWithDefaults(const std::string &name, double value) {
  return getDefaultTables().normalize(name, value);
}
}  // namespace

BOOST_PYTHON_MODULE(rdNormalizedDescriptors) {
  python::scope().attr("__doc__") =
      "Module containing histogram normalization of molecular descriptors, a "
      "port of the RDKit2DHistogramNormalized descriptors from "
      "descriptastorus";

  python::class_<CumulativeHistogram>(
      "CumulativeHistogram",
      "A cumulative histogram for one descriptor: the left edge of each bin "
      "and the fraction of the reference set in that bin or below it.\n"
      "A value normalizes to the fraction of the first bin whose edge is >= "
      "the value, or 1.0 past the last edge. NaN normalizes to 0.0.",
      python::no_init)
      .def("__init__",
           python::make_constructor(
               makeCumulativeHistogram, python::default_call_policies(),
               (python::arg("edges"), python::arg("fractions"))),
           "Constructor. edges must be sorted and the same length as "
           "fractions")
      .def("Normalize", &CumulativeHistogram::normalize,
           (python::arg("self"), python::arg("value")),
           "returns the normalized value")
      .def("__call__", &CumulativeHistogram::normalize,
           (python::arg("self"), python::arg("value")))
      .def("GetEdges", getEdges, python::arg("self"))
      .def("GetFractions", getFractions, python::arg("self"));

  python::class_<HistogramSet>(
      "HistogramSet",
      "A collection of histogram tables keyed by descriptor name",
      python::init<>(python::args("self")))
      .def("AddTable", &HistogramSet::addTable,
           (python::arg("self"), python::arg("name"), python::arg("table")),
           "adds (or replaces) the table for a descriptor")
      .def("HasTable", &HistogramSet::hasTable,
           (python::arg("self"), python::arg("name")))
      .def("GetTable",
           (const CumulativeHistogram &(HistogramSet::*)(size_t) const) &
               HistogramSet::getTable,
           (python::arg("self"), python::arg("idx")),
           python::return_internal_reference<1>(),
           "returns the table with an index, raises IndexError if out of "
           "range")
      .def("GetTable",
           (const CumulativeHistogram &(HistogramSet::*)(const std::string &)
                const) &
               HistogramSet::getTable,
           (python::arg("self"), python::arg("name")),
           python::return_internal_reference<1>(),
           "returns the table for a descriptor, raises KeyError if missing")
      .def("GetNames", getNames, python::arg("self"),
           "returns the descriptor names, in index order")
      .def("GetTableIndex", &HistogramSet::getTableIndex,
           (python::arg("self"), python::arg("name")),
           "returns the index of the table for a descriptor, -1 if missing")
      .def("Normalize",
           (double (HistogramSet::*)(size_t, double) const) &
               HistogramSet::normalize,
           (python::arg("self"), python::arg("idx"), python::arg("value")),
           "normalizes a value with the table with an index")
      .def("Normalize",
           (double (HistogramSet::*)(const std::string &, double) const) &
               HistogramSet::normalize,
           (python::arg("self"), python::arg("name"), python::arg("value")),
           "normalizes a descriptor value, returns 0.0 for descriptors "
           "without a table")
      .def("LoadFromFile", &HistogramSet::loadFromFile,
           (python::arg("self"), python::arg("fileName")),
           "adds the tables in a file to this set")
      .def("LoadFromString", loadFromString,
           (python::arg("self"), python::arg("text")),
           "adds the tables in a string to this set")
      .def("ToString", toString, python::arg("self"),
           "returns the tables in the text format read by LoadFromString")
      .def("__len__", &HistogramSet::size, python::arg("self"));

  // NormalizedProperties derives from rdMolDescriptors.Properties
  python::import("rdkit.Chem.rdMolDescriptors");
  python::class_<NormalizedProperties, boost::shared_ptr<NormalizedProperties>,
                 python::bases<RDKit::Descriptors::Properties>,
                 boost::noncopyable>(
      "NormalizedProperties",
      "Computes descriptors with the rdMolDescriptors.Properties registry and "
      "normalizes them with histogram tables.\n"
      "names defaults to GetNormalizedDescriptorNames() and tables to "
      "GetDefaultTables(). ComputeProperties returns normalized values",
      python::no_init)
      .def("__init__",
           python::make_constructor(makeNormalizedProperties,
                                    python::default_call_policies(),
                                    (python::arg("names") = python::object(),
                                     python::arg("tables") = python::object())))
      .def("GetDescriptorNames", getDescriptorNames, python::arg("self"),
           "returns the descriptor names, in the order values are returned")
      .def("ComputeRawProperties", computeRawHelper,
           (python::arg("self"), python::arg("mol")),
           "returns the unnormalized descriptor values; failures are nan")
      .def("GetTables", &NormalizedProperties::getTables,
           python::return_internal_reference<1>(), python::arg("self"));

  python::def("GetPropertyName", getPropertyName, python::arg("name"),
              "returns the name of the rdMolDescriptors property that computes "
              "a descriptor, e.g. 'exactmw' for 'ExactMolWt'");
  python::def("GetDefaultTablePath", getDefaultTablePath,
              "returns the path of the default histogram tables");
  python::def("GetDefaultTables", getDefaultTables,
              python::return_value_policy<python::reference_existing_object>(),
              "returns the default histogram tables (descriptastorus' "
              "RDKit2DHistogramNormalized reference distributions)");
  python::def("GetNormalizedDescriptorNames",
              getNormalizedDescriptorNamesHelper,
              "returns the names of the descriptors calculated by "
              "CalcNormalizedDescriptors, in the order they are returned "
              "(the order of rdkit.Chem.Descriptors._descList)");
  python::def("CalcNormalizedDescriptors", calcNormalizedHelper,
              (python::arg("mol"), python::arg("tables") = python::object()),
              "calculates the descriptors named by "
              "GetNormalizedDescriptorNames and normalizes them "
              "with tables (the default tables when None). Descriptors that "
              "cannot be calculated or have no table are 0.0");
  python::def("CalcDescriptorValues", calcDescriptorValuesHelper,
              python::arg("mol"),
              "calculates the raw (unnormalized) values of the descriptors "
              "named by GetNormalizedDescriptorNames; failures are nan");
  python::def("NormalizeDescriptor", normalizeWithDefaults,
              (python::arg("name"), python::arg("value")),
              "normalizes a descriptor value with the default tables, returns "
              "0.0 for descriptors without a table");
}
