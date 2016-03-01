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
#include <boost/thread/tss.hpp>
#endif

#include <sstream>

#include "seqs.hpp"
namespace python = boost::python;
using namespace RDKit;

namespace RDKit {
void tossit() { throw IndexErrorException(1); }
}

void rdExceptionTranslator(RDKit::ConformerException const& x) {
  RDUNUSED_PARAM(x);
  PyErr_SetString(PyExc_ValueError, "Bad Conformer Id");
}

void rdSanitExceptionTranslator(RDKit::MolSanitizeException const& x) {
  std::ostringstream ss;
  ss << "Sanitization error: " << x.message();
  PyErr_SetString(PyExc_ValueError, ss.str().c_str());
}

void wrap_table();
void wrap_atom();
void wrap_conformer();
void wrap_bond();
void wrap_mol();
void wrap_ringinfo();
void wrap_EditableMol();
void wrap_monomerinfo();
void wrap_resmolsupplier();

struct PySysErrWrite : std::ostream, std::streambuf
{
  std::string prefix;
  
  PySysErrWrite(const std::string &prefix) :
      std::ostream(this), prefix(prefix) {}

  int overflow(int c) { write(c); return 0;}
  
#ifdef RDK_THREADSAFE_SSS
  void write(char c) { // enable thread safe logging
    static boost::thread_specific_ptr< std::string > buffer;
    if( !buffer.get() ) {
      buffer.reset( new std::string );
    }
    (*buffer.get()) += c;
    if (c == '\n') {
      // Python IO is not thread safe, so grab the GIL
      PyGILState_STATE gstate;
      gstate = PyGILState_Ensure();
      PySys_WriteStderr("%s", (prefix + (*buffer.get())).c_str());
      PyGILState_Release(gstate);
      buffer->clear();
    }
  }
  
#else
  std::string buffer; // unlimited! flushes in endl
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
  static PySysErrWrite debug  ("RDKit DEBUG: ");
  static PySysErrWrite error  ("RDKit ERROR: ");
  static PySysErrWrite info   ("RDKit INFO: ");
  static PySysErrWrite warning("RDKit WARNING: ");
  if( rdDebugLog == NULL || rdInfoLog == NULL || rdErrorLog == NULL || rdWarningLog == NULL ){
    RDLog::InitLogs();
  }
  if( rdDebugLog != NULL ) rdDebugLog->AddTee(debug);
  if( rdInfoLog != NULL ) rdInfoLog->AddTee(info);
  if( rdErrorLog != NULL ) rdErrorLog->AddTee(error);
  if( rdWarningLog != NULL ) rdWarningLog->AddTee(warning);
}

BOOST_PYTHON_MODULE(rdchem) {
  python::scope().attr("__doc__") =
      "Module containing the core chemistry functionality of the RDKit";
  RegisterListConverter<RDKit::Atom*>();
  RegisterListConverter<RDKit::Bond*>();
  rdkit_import_array();
  python::register_exception_translator<RDKit::MolSanitizeException>(
      &rdSanitExceptionTranslator);

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
               1, python::with_custodian_and_ward_postcall<0, 1> >())
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
               1, python::with_custodian_and_ward_postcall<0, 1> >())
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
  wrap_mol();
  wrap_EditableMol();
  wrap_ringinfo();
  wrap_monomerinfo();
  wrap_resmolsupplier();

  //*********************************************
  //
  //  Functions
  //
  //*********************************************

  std::string docString;

  python::def("tossit", tossit);
}
