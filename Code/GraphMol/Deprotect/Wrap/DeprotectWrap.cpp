//
//  Copyright (C) 2020 Brian P Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Deprotect/Deprotect.h>

namespace python = boost::python;
namespace RDKit {
// note: boost::python doesn't support unique_ptr so we convert to shared for
//  python.
boost::shared_ptr<ROMol> DeprotectVectWrap(const ROMol &mol,
                                           const python::object &iterable) {
  std::vector<Deprotect::DeprotectData> deprotections;
  pythonObjectToVect<Deprotect::DeprotectData>(iterable, deprotections);
  auto res = Deprotect::deprotect(mol, deprotections);
  auto m = boost::shared_ptr<ROMol>(res.get());
  res.release();
  return m;
}

boost::shared_ptr<ROMol> DeprotectWrap(const ROMol &mol) {
  auto res = Deprotect::deprotect(mol, Deprotect::getDeprotections());
  auto m = boost::shared_ptr<ROMol>(res.get());
  res.release();
  return m;
}

//! Make a copy so we don't try and change a const vector
std::vector<Deprotect::DeprotectData> GetDeprotectionsWrap() { return Deprotect::getDeprotections(); }
}  // namespace RDKit

struct deprotect_wrap {
  static void wrap() {
    const char *constructor_doc =
        "Construct a new DeprotectData instance.\n"
        "  >>> reaction_class = \"amine\"\n"
        "  >>> reaction_smarts = "
        "\"[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]\"\n"
        "  >>> abbreviation = \"Boc\"\n"
        "  >>> full_name = \"tert-butyloxycarbonyl\"\n"
        "  >>> data = DeprotectData(reaction_class, reaction_smarts, "
        "abbreviation, full_name)\n"
        "  >>> assert data.isValid()\n"
        "\n";

    const char *deprotect_doc_string =
        "DeprotectData class, contains a single deprotection reaction and "
        "information\n"
        "\n"
        " deprotectdata.deprotection_class - functional group being protected\n"
        " deprotectdata.reaction_smarts - reaction smarts used for "
        "deprotection\n"
        " deprotectdata.abbreviation - common abbreviation for the protecting "
        "group\n"
        " deprotectdata.full_name - full name for the protecting group\n"
        "\n"
        "\n";

    python::class_<std::vector<RDKit::Deprotect::DeprotectData>>("DeprotectDataVect")
        .def(
	     python::vector_indexing_suite<std::vector<RDKit::Deprotect::DeprotectData>>());

    python::class_<RDKit::Deprotect::DeprotectData>(
        "DeprotectData", deprotect_doc_string,
        python::init<std::string, std::string, std::string, std::string>(
            constructor_doc))
        .def_readonly("deprotection_class",
                      &RDKit::Deprotect::DeprotectData::deprotection_class)
        .def_readonly("full_name", &RDKit::Deprotect::DeprotectData::full_name)
        .def_readonly("abbreviation", &RDKit::Deprotect::DeprotectData::abbreviation)
        .def_readonly("reaction_smarts", &RDKit::Deprotect::DeprotectData::reaction_smarts)
        .def_readonly("example", &RDKit::Deprotect::DeprotectData::example)
        .def("isValid", &RDKit::Deprotect::DeprotectData::isValid,
             "Returns True if the DeprotectData has a valid reaction");

    python::def("GetDeprotections", &RDKit::GetDeprotectionsWrap,
                "Return the default list of deprotections");

    python::def("Deprotect", &RDKit::DeprotectWrap, python::arg("mol"),
                "Return the deprotected version of the molecule.");

    python::def("Deprotect", &RDKit::DeprotectVectWrap, python::arg("mol"),
                python::arg("deprotections"),
                "Given a list of deprotections, return the deprotected version "
                "of the molecule.");
  }
};

void wrap_deprotect() { deprotect_wrap::wrap(); }
