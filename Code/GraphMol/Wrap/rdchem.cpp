// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdchem_array_API
#include "rdchem.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SanitException.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/python/numpy.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <sstream>
#include <utility>

#include "seqs.hpp"
namespace python = boost::python;
using namespace RDKit;

namespace RDKit {
void tossit() { throw IndexErrorException(1); }
}  // namespace RDKit

void rdExceptionTranslator(RDKit::ConformerException const &) {
  PyErr_SetString(PyExc_ValueError, "Bad Conformer Id");
}

void wrap_table();
void wrap_atom();
void wrap_conformer();
void wrap_bond();
void wrap_stereogroup();
void wrap_mol();
void wrap_ringinfo();
void wrap_EditableMol();
void wrap_monomerinfo();
void wrap_resmolsupplier();
void wrap_molbundle();
void wrap_sgroup();
void wrap_chirality();

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

BOOST_PYTHON_MODULE(rdchem) {
  python::scope().attr("__doc__") =
      "Module containing the core chemistry functionality of the RDKit";
  const bool register_scalar_converters=false;
  boost::python::numpy::initialize(register_scalar_converters);
  RegisterListConverter<RDKit::Atom *>();
  RegisterListConverter<RDKit::Bond *>();
  RegisterListConverter<RDKit::CONFORMER_SPTR>();
  rdkit_import_array();

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

  //*********************************************
  //
  //  Classes
  //
  //*********************************************
  wrap_table();
  wrap_atom();
  wrap_conformer();
  wrap_bond();
  wrap_stereogroup();
  wrap_mol();
  wrap_EditableMol();
  wrap_ringinfo();
  wrap_monomerinfo();
  wrap_resmolsupplier();
  wrap_molbundle();
  wrap_sgroup();
  wrap_chirality();

  //*********************************************
  //
  //  Functions
  //
  //*********************************************

  std::string docString;

  python::def("tossit", tossit);
}
