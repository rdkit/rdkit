//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
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
#include <RDBoost/Wrap_nb.h>

#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
std::string serializeFragCatalog(const RDKit::FragCatalog &self) {
  return self.Serialize();
}
}  // namespace

namespace RDKit {

namespace {

unsigned int GetBitEntryId(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  return self->getIdOfEntryWithBitId(idx);
}

unsigned int GetEntryBitId(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithIdx(idx)->getBitId();
}

std::string GetEntryDescription(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithIdx(idx)->getDescription();
}

std::string GetBitDescription(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithBitId(idx)->getDescription();
}

unsigned int GetEntryOrder(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithIdx(idx)->getOrder();
}

unsigned int GetBitOrder(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithBitId(idx)->getOrder();
}

INT_VECT GetEntryFuncGroupIds(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
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
    throw nb::index_error("index out of range");
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
    throw nb::index_error("index out of range");
  }
  return self->getDownEntryList(idx);
}

DOUBLE_VECT GetBitDiscrims(const FragCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  DOUBLE_VECT res;
  const FragCatalogEntry *entry = self->getEntryWithBitId(idx);
  Subgraphs::DiscrimTuple tmp = entry->getDiscrims();
  res.push_back(std::get<0>(tmp));
  res.push_back(std::get<1>(tmp));
  res.push_back(std::get<2>(tmp));
  return res;
}

}  // namespace

}  // namespace RDKit

void wrap_fragcat(nb::module_ &m) {
  // FIX: none of the functions giving access to the entries in the catalog
  // are being exposed to python
  // right now, adding entries for example should happen through the
  // FragCatGenerator

  nb::class_<RDKit::FragCatalog>(m, "FragCatalog")
      .def(nb::new_([]() { return new RDKit::FragCatalog(); }))
      .def(nb::new_([](RDKit::FragCatParams *fc) {
             return new RDKit::FragCatalog(fc);
           }),
           "params"_a)
      .def(nb::new_([](const nb::bytes pkl) {
             return new RDKit::FragCatalog(std::string(
                 static_cast<const char *>(pkl.data()), pkl.size()));
           }),
           "pickle"_a)
      .def("GetNumEntries", &RDKit::FragCatalog::getNumEntries)
      .def("GetFPLength", &RDKit::FragCatalog::getFPLength)
      .def("GetCatalogParams",
           (RDKit::FragCatParams * (RDKit::FragCatalog::*)()) &
               RDKit::FragCatalog::getCatalogParams,
           nb::rv_policy::reference_internal)
      .def("Serialize",
           [](const RDKit::FragCatalog &self) {
             const auto pkl = self.Serialize();
             return nb::bytes(pkl.c_str(), pkl.size());
           })

      .def("GetBitDescription", &RDKit::GetBitDescription, "idx"_a)
      .def("GetBitOrder", &RDKit::GetBitOrder, "idx"_a)
      .def("GetBitFuncGroupIds", &RDKit::GetBitFuncGroupIds, "idx"_a)
      .def("GetBitEntryId", &RDKit::GetBitEntryId, "idx"_a)

      .def("GetEntryBitId", &RDKit::GetEntryBitId, "idx"_a)
      .def("GetEntryDescription", &RDKit::GetEntryDescription, "idx"_a)
      .def("GetEntryOrder", &RDKit::GetEntryOrder, "idx"_a)
      .def("GetEntryFuncGroupIds", &RDKit::GetEntryFuncGroupIds, "idx"_a)
      .def("GetEntryDownIds", &RDKit::GetEntryDownIds, "idx"_a)

      .def("GetBitDiscrims", &RDKit::GetBitDiscrims, "idx"_a)

      // enable pickle support
      .def("__getstate__",
           getObjectState<RDKit::FragCatalog, serializeFragCatalog>)
      .def("__setstate__", setObjectState<RDKit::FragCatalog>)

      ;
}
