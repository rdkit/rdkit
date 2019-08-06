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
#include <RDBoost/Wrap.h>
#include "rdchem.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SanitException.h>
#include <RDBoost/import_array.h>
#include <RDBoost/iterator_next.h>

#ifdef RDK_THREADSAFE_SSS
// Thread local storage for output buffer for RDKit Logging
#include <thread>
#endif

#include <sstream>
#include <utility>

#include "seqs.hpp"
namespace python = boost::python;
using namespace RDKit;

namespace RDKit {
void tossit() { throw IndexErrorException(1); }
}  // namespace RDKit

void rdExceptionTranslator(RDKit::ConformerException const &x) {
  RDUNUSED_PARAM(x);
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

struct PySysErrWrite : std::ostream, std::streambuf {
  std::string prefix;

  PySysErrWrite(std::string prefix)
      : std::ostream(this), prefix(std::move(prefix)) {}

  int overflow(int c) override {
    write(c);
    return 0;
  }

#ifdef RDK_THREADSAFE_SSS
  void write(char c) {  // enable thread safe logging
    static thread_local std::string buffer = "";
    buffer += c;
    if (c == '\n') {
      // Python IO is not thread safe, so grab the GIL
      PyGILState_STATE gstate;
      gstate = PyGILState_Ensure();
      PySys_WriteStderr("%s", (prefix + buffer).c_str());
      PyGILState_Release(gstate);
      buffer.clear();
    }
  }
#else
  std::string buffer;  // unlimited! flushes in endl
  void write(char c) {
    buffer += c;
    if (c == '\n') {
      PySys_WriteStderr("%s", (prefix + buffer).c_str());
      buffer.clear();
    }
  }
#endif
};

void RDLogError(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdErrorLog) << msg.c_str() << std::endl;
}

void RDLogWarning(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdWarningLog) << msg.c_str() << std::endl;
}

void WrapLogs() {
  static PySysErrWrite debug("RDKit DEBUG: ");
  static PySysErrWrite error("RDKit ERROR: ");
  static PySysErrWrite info("RDKit INFO: ");
  static PySysErrWrite warning("RDKit WARNING: ");
  if (!rdDebugLog || !rdInfoLog || !rdErrorLog || !rdWarningLog) {
    RDLog::InitLogs();
  }
  if (rdDebugLog != nullptr) rdDebugLog->SetTee(debug);
  if (rdInfoLog != nullptr) rdInfoLog->SetTee(info);
  if (rdErrorLog != nullptr) rdErrorLog->SetTee(error);
  if (rdWarningLog != nullptr) rdWarningLog->SetTee(warning);
}

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
  PyErr_SetString(pyExcType, x.message());
}

// pattern from here:
// https://stackoverflow.com/questions/9620268/boost-python-custom-exception-class
PyObject *createExceptionClass(const char *name,
                               PyObject *baseTypeObj = PyExc_ValueError) {
  std::string scopeName =
      python::extract<std::string>(python::scope().attr("__name__"));
  std::string qualifiedName0 = scopeName + "." + name;
  char *qualifiedName1 = const_cast<char *>(qualifiedName0.c_str());

  PyObject *typeObj = PyErr_NewException(qualifiedName1, baseTypeObj, 0);
  if (!typeObj) python::throw_error_already_set();
  python::scope().attr(name) = python::handle<>(python::borrowed(typeObj));
  return typeObj;
}

BOOST_PYTHON_MODULE(rdchem) {
  python::scope().attr("__doc__") =
      "Module containing the core chemistry functionality of the RDKit";
  RegisterListConverter<RDKit::Atom *>();
  RegisterListConverter<RDKit::Bond *>();
  rdkit_import_array();

  // this is one of those parts where I think I wish that I knew how to do
  // template meta-programming
  python::class_<MolSanitizeException>("_cppMolSanitizeException",
                                       "exception arising from sanitization",
                                       python::no_init)
      .def("Message", &MolSanitizeException::message)
      .def("GetType", &MolSanitizeException::getType);
  python::register_ptr_to_python<boost::shared_ptr<MolSanitizeException>>();
  molSanitizeExceptionType = createExceptionClass("MolSanitizeException");
  python::register_exception_translator<RDKit::MolSanitizeException>(
      [&](const MolSanitizeException &exc) {
        sanitExceptionTranslator(exc, molSanitizeExceptionType);
      });

  python::class_<AtomSanitizeException, python::bases<MolSanitizeException>>(
      "_cppAtomSanitizeException", "exception arising from sanitization",
      python::no_init)
      .def("GetAtomIdx", &AtomSanitizeException::getAtomIdx);
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

  python::def("WrapLogs", WrapLogs,
              "Wrap the internal RDKit streams so they go to python's "
              "SysStdErr");
  python::def("LogWarningMsg", RDLogWarning,
              "Log a warning message to the RDKit warning logs");
  python::def("LogErrorMsg", RDLogError,
              "Log a warning message to the RDKit error logs");

  //*********************************************
  //
  //  Utility Classes
  //
  //*********************************************
  python::class_<AtomIterSeq>(
      "_ROAtomSeq",
      "Read-only sequence of atoms, not constructable from Python.",
      python::no_init)
      .def("__iter__", &AtomIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def(NEXT_METHOD, &AtomIterSeq::next,
           python::return_value_policy<python::reference_existing_object>())

      .def("__len__", &AtomIterSeq::len)
      .def("__getitem__", &AtomIterSeq::get_item,
           python::return_value_policy<python::reference_existing_object>());
  python::class_<QueryAtomIterSeq>("_ROQAtomSeq",
                                   "Read-only sequence of atoms matching a "
                                   "query, not constructable from Python.",
                                   python::no_init)
      .def("__iter__", &QueryAtomIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def(NEXT_METHOD, &QueryAtomIterSeq::next,
           python::return_value_policy<python::reference_existing_object>())
      .def("__len__", &QueryAtomIterSeq::len)
      .def("__getitem__", &QueryAtomIterSeq::get_item,
           python::return_value_policy<python::reference_existing_object>());
  python::class_<BondIterSeq>(
      "_ROBondSeq",
      "Read-only sequence of bonds, not constructable from Python.",
      python::no_init)
      // FIX: we ought to be able to expose an iteration interface
      .def("__len__", &BondIterSeq::len)
      .def("__getitem__", &BondIterSeq::get_item,
           python::return_value_policy<python::reference_existing_object>());

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

  //*********************************************
  //
  //  Functions
  //
  //*********************************************

  std::string docString;

  python::def("tossit", tossit);
}
