//
//  Copyright (C) 2006-2026 Greg Landrum and other RDKit contributors
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
      .def(
          "__init__",
          [](MolCatalog &self, nb::bytes pickle) {
            new (&self) MolCatalog(std::string(
                static_cast<const char *>(pickle.data()), pickle.size()));
          },
          "pickle"_a)
      .def("GetNumEntries", &MolCatalog::getNumEntries)
      .def("GetFPLength", &MolCatalog::getFPLength)
      .def("Serialize",
           [](const MolCatalog &self) {
             const auto pkl = self.Serialize();
             return nb::bytes(pkl.c_str(), pkl.size());
           })
      .def("GetBitDescription", GetBitDescription, "idx"_a)
      .def("GetBitEntryId", GetBitEntryId, "idx"_a)
      .def("GetEntryBitId", GetEntryBitId, "idx"_a)
      .def("GetEntryDescription", GetEntryDescription, "idx"_a)
      .def("GetEntryDownIds", GetEntryDownIds, "idx"_a)
      .def("AddEntry", AddEntry, "entry"_a)
      .def("AddEdge", &MolCatalog::addEdge, "id1"_a, "id2"_a)
      // pickle support
      .def("__getstate__",
           [](const MolCatalog &self) {
             const auto pkl = self.Serialize();
             return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalog &self, const std::tuple<nb::bytes> &state) {
             const auto &pkl = std::get<0>(state);
             new (&self) MolCatalog(std::string(
                 static_cast<const char *>(pkl.data()), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalog &self, const std::tuple<std::string> &state) {
             new (&self) MolCatalog(std::get<0>(state));
           });

  nb::class_<MolCatalogEntry>(m, "MolCatalogEntry")
      .def(nb::init<>())
      .def(
          "__init__",
          [](MolCatalogEntry &self, nb::bytes pickle) {
            new (&self) MolCatalogEntry(std::string(
                static_cast<const char *>(pickle.data()), pickle.size()));
          },
          "pickle"_a)
      .def("GetDescription", &MolCatalogEntry::getDescription)
      .def("SetDescription", &MolCatalogEntry::setDescription, "val"_a)
      .def("GetMol", catalogEntryGetMol, nb::rv_policy::reference_internal)
      .def("SetMol", catalogEntrySetMol, "mol"_a)
      .def("GetOrder", &MolCatalogEntry::getOrder)
      .def("SetOrder", &MolCatalogEntry::setOrder, "order"_a)
      .def("Serialize",
           [](const MolCatalogEntry &self) {
             const auto pkl = self.Serialize();
             return nb::bytes(pkl.c_str(), pkl.size());
           })
      // pickle support
      .def("__getstate__",
           [](const MolCatalogEntry &self) {
             const auto pkl = self.Serialize();
             return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalogEntry &self, const std::tuple<nb::bytes> &state) {
             const auto &pkl = std::get<0>(state);
             new (&self) MolCatalogEntry(std::string(
                 static_cast<const char *>(pkl.data()), pkl.size()));
           })
      .def("__setstate__",
           [](MolCatalogEntry &self, const std::tuple<std::string> &state) {
             new (&self) MolCatalogEntry(std::get<0>(state));
           });

  m.def("CreateMolCatalog", createMolCatalog, nb::rv_policy::take_ownership);
}
