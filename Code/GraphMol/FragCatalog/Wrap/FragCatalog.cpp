//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>

namespace python = boost::python;
namespace RDKit {

struct fragcatalog_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const FragCatalog &self) {
    std::string res;
    res = self.Serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};
unsigned int GetBitEntryId(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  return self->getIdOfEntryWithBitId(idx);
}

unsigned int GetEntryBitId(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getEntryWithIdx(idx)->getBitId();
}
std::string GetEntryDescription(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getEntryWithIdx(idx)->getDescription();
}
std::string GetBitDescription(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  return self->getEntryWithBitId(idx)->getDescription();
}
unsigned int GetEntryOrder(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getEntryWithIdx(idx)->getOrder();
}
unsigned int GetBitOrder(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  return self->getEntryWithBitId(idx)->getOrder();
}
INT_VECT GetEntryFuncGroupIds(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  INT_VECT res;
  INT_INT_VECT_MAP gps = self->getEntryWithIdx(idx)->getFuncGroupMap();
  for (const auto &iv : gps) {
    for (int ivci : iv.second) {
      res.push_back(ivci);
    }
  }
  return res;
}
INT_VECT GetBitFuncGroupIds(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  INT_VECT res;
  INT_INT_VECT_MAP gps = self->getEntryWithBitId(idx)->getFuncGroupMap();
  for (const auto &iv : gps) {
    for (int ivci : iv.second) {
      res.push_back(ivci);
    }
  }
  return res;
}
INT_VECT GetEntryDownIds(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getDownEntryList(idx);
}

DOUBLE_VECT GetBitDiscrims(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  DOUBLE_VECT res;
  const FragCatalogEntry *entry = self->getEntryWithBitId(idx);
  Subgraphs::DiscrimTuple tmp = entry->getDiscrims();
  res.push_back(std::get<0>(tmp));
  res.push_back(std::get<1>(tmp));
  res.push_back(std::get<2>(tmp));
  return res;
}

struct fragcat_wrapper {
  static void wrap() {
    // FIX: none of the functions giving access to the entries in the catalog
    // are  being exposed to python
    // right now, adding entries for example should happen through the
    // FragCatGenerator

    python::class_<FragCatalog>(
        "FragCatalog",
        python::init<FragCatParams *>(python::args("self", "params")))
        .def(python::init<const std::string &>(python::args("self", "pickle")))
        .def("GetNumEntries", &FragCatalog::getNumEntries, python::args("self"))
        .def("GetFPLength", &FragCatalog::getFPLength, python::args("self"))
        .def("GetCatalogParams",
             (FragCatParams * (FragCatalog::*)()) &
                 FragCatalog::getCatalogParams,
             python::return_value_policy<python::reference_existing_object>(),
             python::args("self"))
        .def("Serialize", &FragCatalog::Serialize, python::args("self"))

        .def("GetBitDescription", &GetBitDescription,
             python::args("self", "idx"))
        .def("GetBitOrder", &GetBitOrder, python::args("self", "idx"))
        .def("GetBitFuncGroupIds", &GetBitFuncGroupIds,
             python::args("self", "idx"))
        .def("GetBitEntryId", &GetBitEntryId, python::args("self", "idx"))

        .def("GetEntryBitId", &GetEntryBitId, python::args("self", "idx"))
        .def("GetEntryDescription", &GetEntryDescription,
             python::args("self", "idx"))
        .def("GetEntryOrder", &GetEntryOrder, python::args("self", "idx"))
        .def("GetEntryFuncGroupIds", &GetEntryFuncGroupIds,
             python::args("self", "idx"))
        .def("GetEntryDownIds", &GetEntryDownIds, python::args("self", "idx"))

        .def("GetBitDiscrims", &GetBitDiscrims, python::args("self", "idx"))

        // enable pickle support
        .def_pickle(fragcatalog_pickle_suite())

        ;
  };
};

}  // namespace RDKit

void wrap_fragcat() { RDKit::fragcat_wrapper::wrap(); }
