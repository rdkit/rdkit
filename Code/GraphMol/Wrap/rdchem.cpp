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

// std::ostream wrapper around Python's stderr stream
struct PyErrStream : std::ostream, std::streambuf {
  std::string prefix;

  PyErrStream(std::string prefix)
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
      {
        PyGILStateHolder h;
        PySys_WriteStderr("%s%s", prefix.c_str(), buffer.c_str());
      }
      buffer.clear();
    }
  }
#else
  std::string buffer;  // unlimited! flushes in endl
  void write(char c) {
    buffer += c;
    if (c == '\n') {
      PySys_WriteStderr("%s%s", prefix.c_str(), buffer.c_str());
      buffer.clear();
    }
  }
#endif
};

// std::ostream wrapper around Python's logging module
struct PyLogStream : std::ostream, std::streambuf {
  PyObject *logfn = nullptr;
#ifdef RDK_THREADSAFE_SSS
  static thread_local std::string buffer;
#else
  std::string buffer = "";
#endif

  PyLogStream(std::string level): std::ostream(this) {
    PyObject *module = PyImport_ImportModule("logging");
    PyObject *logger = nullptr;

    if (module != nullptr) {
      logger = PyObject_CallMethod(module, "getLogger", "s", "rdkit");
      Py_DECREF(module);
    }

    if (logger != nullptr) {
      logfn = PyObject_GetAttrString(logger, level.c_str());
      Py_DECREF(logger);
    }

    if (PyErr_Occurred()) {
      PyErr_Print();
    }
  }

  ~PyLogStream() {
    if (!_Py_IsFinalizing()) {
      Py_XDECREF(logfn);
    }
  }

  int overflow(int c) override {
    write(c);
    return 0;
  }

  void write(char c) {
    if (logfn == nullptr) {
      return;
    }

    if (c == '\n') {
      PyObject *message = PyUnicode_FromString(buffer.c_str());
      if (message != nullptr) {
#ifdef RDK_THREADSAFE_SSS
        PyGILStateHolder h;
#endif
        PyObject *result = PyObject_CallOneArg(logfn, message);
        Py_XDECREF(result);
        Py_DECREF(message);
      }

      buffer.clear();
    }
    else {
      buffer += c;
    }
  }
};

#ifdef RDK_THREADSAFE_SSS
thread_local std::string PyLogStream::buffer;
#endif

void RDLogDebug(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdDebugLog) << msg.c_str() << std::endl;
}

void RDLogInfo(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdInfoLog) << msg.c_str() << std::endl;
}

void RDLogWarning(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdWarningLog) << msg.c_str() << std::endl;
}

void RDLogError(const std::string &msg) {
  NOGIL gil;
  BOOST_LOG(rdErrorLog) << msg.c_str() << std::endl;
}

void WrapLogs() {
  static PyErrStream debug("RDKit DEBUG: ");
  static PyErrStream error("RDKit ERROR: ");
  static PyErrStream info("RDKit INFO: ");
  static PyErrStream warning("RDKit WARNING: ");
  if (!rdDebugLog || !rdInfoLog || !rdErrorLog || !rdWarningLog) {
    RDLog::InitLogs();
  }
  if (rdDebugLog != nullptr) {
    rdDebugLog->SetTee(debug);
  }
  if (rdInfoLog != nullptr) {
    rdInfoLog->SetTee(info);
  }
  if (rdErrorLog != nullptr) {
    rdErrorLog->SetTee(error);
  }
  if (rdWarningLog != nullptr) {
    rdWarningLog->SetTee(warning);
  }
}

void LogToPythonLogger() {
  static PyLogStream debug("debug");
  static PyLogStream info("info");
  static PyLogStream warning("warning");
  static PyLogStream error("error");

  rdDebugLog   = std::make_shared<boost::logging::rdLogger>(&debug);
  rdInfoLog    = std::make_shared<boost::logging::rdLogger>(&info);
  rdWarningLog = std::make_shared<boost::logging::rdLogger>(&warning);
  rdErrorLog   = std::make_shared<boost::logging::rdLogger>(&error);
}

void LogToPythonStderr() {
  static PyErrStream debug("");
  static PyErrStream info("");
  static PyErrStream warning("");
  static PyErrStream error("");

  rdDebugLog   = std::make_shared<boost::logging::rdLogger>(&debug);
  rdInfoLog    = std::make_shared<boost::logging::rdLogger>(&info);
  rdWarningLog = std::make_shared<boost::logging::rdLogger>(&warning);
  rdErrorLog   = std::make_shared<boost::logging::rdLogger>(&error);
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
  RegisterListConverter<RDKit::Atom *>();
  RegisterListConverter<RDKit::Bond *>();
  RegisterListConverter<RDKit::CONFORMER_SPTR>();
  rdkit_import_array();

  // this is one of those parts where I think I wish that I knew how to do
  // template meta-programming
  python::class_<MolSanitizeException>("_cppMolSanitizeException",
                                       "exception arising from sanitization",
                                       python::no_init)
      .def("Message", &MolSanitizeException::what)
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

  python::def("LogToCppStreams", RDLog::InitLogs,
              "Initialize RDKit logs with C++ streams");
  python::def("LogToPythonLogger", LogToPythonLogger,
              "Initialize RDKit logs with Python's logging module");
  python::def("LogToPythonStderr", LogToPythonStderr,
              "Initialize RDKit logs with Python's stderr stream");
  python::def("WrapLogs", WrapLogs,
              "Tee RDKit logs to Python's stderr stream");

  python::def("LogDebugMsg", RDLogDebug,
              "Log a message to the RDKit debug logs");
  python::def("LogInfoMsg", RDLogInfo,
              "Log a message to the RDKit info logs");
  python::def("LogWarningMsg", RDLogWarning,
              "Log a message to the RDKit warning logs");
  python::def("LogErrorMsg", RDLogError,
              "Log a message to the RDKit error logs");

  //*********************************************
  //
  //  Utility Classes
  //
  //*********************************************
  python::class_<AtomIterSeq>(
      "_ROAtomSeq",
      "Read-only sequence of atoms, not constructible from Python.",
      python::no_init)
      .def("__iter__", &AtomIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__next__", &AtomIterSeq::next,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())

      .def("__len__", &AtomIterSeq::len)
      .def("__getitem__", &AtomIterSeq::get_item,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>());
  python::class_<QueryAtomIterSeq>("_ROQAtomSeq",
                                   "Read-only sequence of atoms matching a "
                                   "query, not constructible from Python.",
                                   python::no_init)
      .def("__iter__", &QueryAtomIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__next__", &QueryAtomIterSeq::next,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__len__", &QueryAtomIterSeq::len)
      .def("__getitem__", &QueryAtomIterSeq::get_item,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>());
  python::class_<BondIterSeq>(
      "_ROBondSeq",
      "Read-only sequence of bonds, not constructible from Python.",
      python::no_init)
      .def("__iter__", &BondIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__next__", &BondIterSeq::next,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__len__", &BondIterSeq::len)
      .def("__getitem__", &BondIterSeq::get_item,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>());
  python::class_<ConformerIterSeq>(
      "_ROConformerSeq",
      "Read-only sequence of conformers, not constructible from Python.",
      python::no_init)
      .def("__iter__", &ConformerIterSeq::__iter__,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())
      .def("__next__", next_ptr<ConformerIterSeq, Conformer>,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>())

      .def("__len__", &ConformerIterSeq::len)
      .def("__getitem__", get_item_ptr<ConformerIterSeq, Conformer>,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>());

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
