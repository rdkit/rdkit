//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include <sstream>
#include <utility>

// #include "seqs.hpp"
using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

void wrap_table(nb::module_ &);
void wrap_mol(nb::module_ &);
void wrap_atom(nb::module_ &);
void wrap_bond(nb::module_ &);
// void wrap_conformer(nb::module_ &);
// void wrap_stereogroup(nb::module_ &);
// void wrap_ringinfo(nb::module_ &);
// void wrap_EditableMol(nb::module_ &);
// void wrap_monomerinfo(nb::module_ &);
// void wrap_resmolsupplier(nb::module_ &);
// void wrap_molbundle(nb::module_ &);
// void wrap_sgroup(nb::module_ &);
// void wrap_chirality(nb::module_ &);
#if 0
python::tuple getAtomIndicesHelper(const KekulizeException &self) {
  python::list res;
  for (auto idx : self.getAtomIndices()) {
    res.append(idx);
  }
  return python::tuple(res);
}

PyObject *molSanitizeExceptionType = nullptr;
PyObject *atomSanitizeExceptionType = nullptr;
PyObject *atomValenceExceptionType = nullptr;
PyObject *atomKekulizeExceptionType = nullptr;
PyObject *kekulizeExceptionType = nullptr;

// pattern from here:
// https://stackoverflow.com/questions/11448735/boostpython-export-custom-exception-and-inherit-from-pythons-exception
template <typename EXC_TYPE>
void sanitExceptionTranslator(const EXC_TYPE &x, PyObject *pyExcType) {
  PRECONDITION(pyExcType != nullptr, "global type not initialized");
  python::object pyExcInstance(python::handle<>(python::borrowed(pyExcType)));
  pyExcInstance.attr("cause") = x;
  PyErr_SetString(pyExcType, x.what());
}

// pattern from here:
// https://stackoverflow.com/questions/9620268/boost-python-custom-exception-class
PyObject *createExceptionClass(const char *name,
                               PyObject *baseTypeObj = PyExc_ValueError) {
  std::string scopeName =
      python::extract<std::string>(python::scope().attr("__name__"));
  std::string qualifiedName0 = scopeName + "." + name;
  char *qualifiedName1 = const_cast<char *>(qualifiedName0.c_str());

  PyObject *typeObj = PyErr_NewException(qualifiedName1, baseTypeObj, nullptr);
  if (!typeObj) {
    python::throw_error_already_set();
  }
  python::scope().attr(name) = python::handle<>(python::borrowed(typeObj));
  return typeObj;
}

template <typename O, typename T>
T *get_item_ptr(O &self, int i) {
  return self.get_item(i).get();
}

template <typename O, typename T>
T *next_ptr(O &self) {
  return self.next().get();
}
#endif
NB_MODULE(rdchem, m) {
  m.doc() = "Module containing the core chemistry functionality of the RDKit";
  // const bool register_scalar_converters = false;
  // boost::python::numpy::initialize(register_scalar_converters);
  // RegisterListConverter<RDKit::Atom *>();
  // RegisterListConverter<RDKit::Bond *>();
  // RegisterListConverter<RDKit::CONFORMER_SPTR>();
  // rdkit_import_array();

#if 0
  // this is one of those parts where I think I wish that I knew how to do
  // template meta-programming
  python::class_<MolSanitizeException>("_cppMolSanitizeException",
                                       "exception arising from sanitization",
                                       python::no_init)
      .def("Message", &MolSanitizeException::what, python::args("self"))
      .def("GetType", &MolSanitizeException::getType, python::args("self"));
  python::register_ptr_to_python<boost::shared_ptr<MolSanitizeException>>();
  molSanitizeExceptionType = createExceptionClass("MolSanitizeException");
  python::register_exception_translator<RDKit::MolSanitizeException>(
      [&](const MolSanitizeException &exc) {
        sanitExceptionTranslator(exc, molSanitizeExceptionType);
      });

  python::class_<AtomSanitizeException, python::bases<MolSanitizeException>>(
      "_cppAtomSanitizeException", "exception arising from sanitization",
      python::no_init)
      .def("GetAtomIdx", &AtomSanitizeException::getAtomIdx,
           python::args("self"));
  python::register_ptr_to_python<boost::shared_ptr<AtomSanitizeException>>();
  atomSanitizeExceptionType =
      createExceptionClass("AtomSanitizeException", molSanitizeExceptionType);
  python::register_exception_translator<RDKit::AtomSanitizeException>(
      [&](const AtomSanitizeException &exc) {
        sanitExceptionTranslator(exc, atomSanitizeExceptionType);
      });

  python::class_<AtomValenceException, python::bases<AtomSanitizeException>>(
      "_cppAtomValenceException", "exception arising from sanitization",
      python::no_init);
  python::register_ptr_to_python<boost::shared_ptr<AtomValenceException>>();
  atomValenceExceptionType =
      createExceptionClass("AtomValenceException", atomSanitizeExceptionType);
  python::register_exception_translator<RDKit::AtomValenceException>(
      [&](const AtomValenceException &exc) {
        sanitExceptionTranslator(exc, atomValenceExceptionType);
      });

  python::class_<AtomKekulizeException, python::bases<AtomSanitizeException>>(
      "_cppAtomKekulizeException", "exception arising from sanitization",
      python::no_init);
  python::register_ptr_to_python<boost::shared_ptr<AtomKekulizeException>>();
  atomKekulizeExceptionType =
      createExceptionClass("AtomKekulizeException", atomSanitizeExceptionType);
  python::register_exception_translator<RDKit::AtomKekulizeException>(
      [&](const AtomKekulizeException &exc) {
        sanitExceptionTranslator(exc, atomKekulizeExceptionType);
      });

  python::class_<KekulizeException, python::bases<MolSanitizeException>>(
      "_cppAtomKekulizeException", "exception arising from sanitization",
      python::no_init)
      .def("GetAtomIndices", &getAtomIndicesHelper);
  python::register_ptr_to_python<boost::shared_ptr<KekulizeException>>();
  kekulizeExceptionType =
      createExceptionClass("KekulizeException", molSanitizeExceptionType);
  python::register_exception_translator<RDKit::KekulizeException>(
      [&](const KekulizeException &exc) {
        sanitExceptionTranslator(exc, kekulizeExceptionType);
      });
      //*********************************************
      //
      //  Utility Classes
      //
      //*********************************************
      python::class_<QueryAtomIterSeq>("_ROQAtomSeq",
        "Read-only sequence of atoms matching a "
        "query, not constructible from Python.",
        python::no_init)
        .def("__iter__", &QueryAtomIterSeq::__iter__,
          python::return_internal_reference<
          1, python::with_custodian_and_ward_postcall<0, 1>>(),
          python::args("self"))
          .def("__next__", &QueryAtomIterSeq::next,
            python::return_internal_reference<
            1, python::with_custodian_and_ward_postcall<0, 1>>(),
            python::args("self"))
            .def("__len__", &QueryAtomIterSeq::len, python::args("self"))
            .def("__getitem__", &QueryAtomIterSeq::get_item,
              python::return_internal_reference<
              1, python::with_custodian_and_ward_postcall<0, 1>>(),
              python::args("self", "which"));
              python::class_<ConformerIterSeq>(
                "_ROConformerSeq",
                "Read-only sequence of conformers, not constructible from Python.",
                python::no_init)
                .def("__iter__", &ConformerIterSeq::__iter__,
                  python::return_internal_reference<
                  1, python::with_custodian_and_ward_postcall<0, 1>>(),
                  python::args("self"))
                  .def("__next__", next_ptr<ConformerIterSeq, Conformer>,
                    python::return_internal_reference<
                    1, python::with_custodian_and_ward_postcall<0, 1>>(),
                    python::args("self"))
                    
                    .def("__len__", &ConformerIterSeq::len, python::args("self"))
                    .def("__getitem__", get_item_ptr<ConformerIterSeq, Conformer>,
                      python::return_internal_reference<
                      1, python::with_custodian_and_ward_postcall<0, 1>>(),
                      python::args("self", "i"));
#endif

  //*********************************************
  //
  //  Classes
  //
  //*********************************************
  wrap_table(m);
  wrap_mol(m);
  wrap_atom(m);
  // wrap_bond(m);
  // wrap_stereogroup();
  // wrap_conformer();
  // wrap_EditableMol();
  // wrap_ringinfo();
  // wrap_monomerinfo();
  // wrap_resmolsupplier();
  // wrap_molbundle();
  // wrap_sgroup();
  // wrap_chirality();

  //*********************************************
  //
  //  Functions
  //
  //*********************************************
  m.def(
      "MolFromSmiles",
      [](const std::string &smiles) {
        auto molp = v2::SmilesParse::MolFromSmiles(smiles);
        std::unique_ptr<ROMol> molp2(molp.release());
        return molp2;
      },
      "smiles"_a, "Constructs a molecule from a SMILES string");
}
