//
//  Copyright (C) 2006 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolCatalog/MolCatalog.h>
#include <GraphMol/MolCatalog/MolCatalogEntry.h>
#include <GraphMol/MolCatalog/MolCatalogParams.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

unsigned int GetBitEntryId(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  return self->getIdOfEntryWithBitId(idx);
}

unsigned int GetEntryBitId(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithIdx(idx)->getBitId();
}

std::string GetEntryDescription(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithIdx(idx)->getDescription();
}

std::string GetBitDescription(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw nb::index_error("index out of range");
  }
  return self->getEntryWithBitId(idx)->getDescription();
}

INT_VECT GetEntryDownIds(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw nb::index_error("index out of range");
  }
  return self->getDownEntryList(idx);
}

unsigned int AddEntry(MolCatalog *self, MolCatalogEntry *entry) {
  auto *cpy = new MolCatalogEntry(*entry);
  return self->addEntry(cpy);
}

void catalogEntrySetMol(MolCatalogEntry *self, const ROMol *mol) {
  auto *cpy = new ROMol(*mol);
  self->setMol(cpy);
}

const ROMol &catalogEntryGetMol(MolCatalogEntry &self) {
  return *self.getMol();
}

MolCatalog *createMolCatalog() {
  MolCatalogParams params;
  return new MolCatalog(&params);
}

}  // namespace

NB_MODULE(rdMolCatalog, m) {
  nb::class_<MolCatalog>(m, "MolCatalog")
      .def(nb::init<const std::string &>(), "pickle"_a)
      .def("GetNumEntries", &MolCatalog::getNumEntries)
      .def("GetFPLength", &MolCatalog::getFPLength)
      .def("Serialize", &MolCatalog::Serialize)
      .def("GetBitDescription", GetBitDescription, "self"_a, "idx"_a)
      .def("GetBitEntryId", GetBitEntryId, "self"_a, "idx"_a)
      .def("GetEntryBitId", GetEntryBitId, "self"_a, "idx"_a)
      .def("GetEntryDescription", GetEntryDescription, "self"_a, "idx"_a)
      .def("GetEntryDownIds", GetEntryDownIds, "self"_a, "idx"_a)
      .def("AddEntry", AddEntry, "self"_a, "entry"_a)
      .def("AddEdge", &MolCatalog::addEdge, "self"_a, "id1"_a, "id2"_a)
      // pickle support
      .def("__getstate__",
           [](const MolCatalog &self) {
             const auto pkl = self.Serialize();
             return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalog &self, const std::tuple<nb::bytes> &state) {
             std::string pkl(
                 static_cast<const char *>(std::get<0>(state).data()),
                 static_cast<size_t>(std::get<0>(state).size()));
             new (&self) MolCatalog(pkl);
           });

  nb::class_<MolCatalogEntry>(m, "MolCatalogEntry")
      .def(nb::init<>())
      .def(nb::init<const std::string &>(), "pickle"_a)
      .def("GetDescription", &MolCatalogEntry::getDescription)
      .def("SetDescription", &MolCatalogEntry::setDescription, "val"_a)
      .def("GetMol", catalogEntryGetMol, nb::rv_policy::reference_internal)
      .def("SetMol", catalogEntrySetMol, "mol"_a)
      .def("GetOrder", &MolCatalogEntry::getOrder)
      .def("SetOrder", &MolCatalogEntry::setOrder, "order"_a)
      // pickle support
      .def("__getstate__",
           [](const MolCatalogEntry &self) {
             const auto pkl = self.Serialize();
             return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalogEntry &self, const std::tuple<nb::bytes> &state) {
             std::string pkl(
                 static_cast<const char *>(std::get<0>(state).data()),
                 static_cast<size_t>(std::get<0>(state).size()));
             new (&self) MolCatalogEntry(pkl);
           });

  m.def("CreateMolCatalog", createMolCatalog, nb::rv_policy::take_ownership);
}
