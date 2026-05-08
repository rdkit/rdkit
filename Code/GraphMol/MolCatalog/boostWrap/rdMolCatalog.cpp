// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
#include "rdMolCatalog.h"
#include <RDBoost/python.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolCatalog/MolCatalog.h>
#include <GraphMol/MolCatalog/MolCatalogEntry.h>
#include <GraphMol/MolCatalog/MolCatalogParams.h>

namespace python = boost::python;
using namespace RDKit;
namespace {
struct molcatalog_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const MolCatalog &self) {
    std::string res;
    res = self.Serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

struct molcatalogentry_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const MolCatalogEntry &self) {
    std::string res;
    res = self.Serialize();
    return python::make_tuple(python::object(python::handle<>(
        PyBytes_FromStringAndSize(res.c_str(), res.length()))));
  };
};

unsigned int GetBitEntryId(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  return self->getIdOfEntryWithBitId(idx);
}

unsigned int GetEntryBitId(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getEntryWithIdx(idx)->getBitId();
}
std::string GetEntryDescription(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getEntryWithIdx(idx)->getDescription();
}
std::string GetBitDescription(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getFPLength()) {
    throw_index_error(idx);
  }
  return self->getEntryWithBitId(idx)->getDescription();
}
INT_VECT GetEntryDownIds(const MolCatalog *self, unsigned int idx) {
  if (idx > self->getNumEntries()) {
    throw_index_error(idx);
  }
  return self->getDownEntryList(idx);
}

unsigned int AddEntry(MolCatalog *self, MolCatalogEntry *entry) {
  auto *cpy = new MolCatalogEntry(*entry);
  return self->addEntry(cpy);
  // return self->addEntry(entry);
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
struct MolCatalog_wrapper {
  static void wrap() {
    python::class_<MolCatalog>(
        "MolCatalog",
        python::init<const std::string &>(python::args("self", "pickle")))
        .def("GetNumEntries", &MolCatalog::getNumEntries, python::args("self"))
        .def("GetFPLength", &MolCatalog::getFPLength, python::args("self"))
        .def("Serialize", &MolCatalog::Serialize, python::args("self"))

        .def("GetBitDescription", GetBitDescription,
             python::args("self", "idx"))
        .def("GetBitEntryId", GetBitEntryId, python::args("self", "idx"))

        .def("GetEntryBitId", GetEntryBitId, python::args("self", "idx"))
        .def("GetEntryDescription", GetEntryDescription,
             python::args("self", "idx"))
        .def("GetEntryDownIds", GetEntryDownIds, python::args("self", "idx"))

        .def("AddEntry", AddEntry, python::args("self", "entry"))
        .def("AddEdge", &MolCatalog::addEdge,
             python::args("self", "id1", "id2"))

        // enable pickle support
        .def_pickle(molcatalog_pickle_suite());
    python::def("CreateMolCatalog", createMolCatalog,
                python::return_value_policy<python::manage_new_object>());
  };
};
struct MolCatalogEntry_wrapper {
  static void wrap() {
    python::class_<MolCatalogEntry>("MolCatalogEntry",
                                    python::init<>(python::args("self")))
        .def(python::init<const std::string &>(python::args("self", "pickle")))
        .def("GetDescription", &MolCatalogEntry::getDescription,
             python::args("self"))
        .def("SetDescription", &MolCatalogEntry::setDescription,
             python::args("self", "val"))
        .def("GetMol", catalogEntryGetMol,
             python::return_internal_reference<1>(), python::args("self"))
        .def("SetMol", catalogEntrySetMol, python::args("self", "mol"))
        .def("GetOrder", &MolCatalogEntry::getOrder, python::args("self"))
        .def("SetOrder", &MolCatalogEntry::setOrder,
             python::args("self", "order"))

        // enable pickle support
        .def_pickle(molcatalogentry_pickle_suite())

        ;
  };
};
}  // namespace

BOOST_PYTHON_MODULE(rdMolCatalog) {
  MolCatalog_wrapper::wrap();
  MolCatalogEntry_wrapper::wrap();
}
