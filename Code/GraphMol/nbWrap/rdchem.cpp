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
#include <nanobind/stl/vector.h>
#include <nanobind/eval.h>

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
void wrap_conformer(nb::module_ &);
void wrap_stereogroup(nb::module_ &);
void wrap_ringinfo(nb::module_ &);
void wrap_EditableMol(nb::module_ &);
void wrap_monomerinfo(nb::module_ &);
void wrap_resmolsupplier(nb::module_ &);
void wrap_molbundle(nb::module_ &);
void wrap_sgroup(nb::module_ &);
void wrap_chirality(nb::module_ &);
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

namespace RDKit {
void tossit() { throw IndexErrorException(1); }
}  // namespace RDKit

NB_MODULE(rdchem, m) {
  m.doc() = "Module containing the core chemistry functionality of the RDKit";
  // const bool register_scalar_converters = false;
  // boost::python::numpy::initialize(register_scalar_converters);
  // RegisterListConverter<RDKit::Atom *>();
  // RegisterListConverter<RDKit::Bond *>();
  // RegisterListConverter<RDKit::CONFORMER_SPTR>();
  // rdkit_import_array();
  m.def("tossit", &RDKit::tossit);
  nb::object scope = m.attr("__dict__");

  nb::exec(R"PY(
class _BaseSanitException(ValueError):
    def __init__(self, msg: str):
        self.cause = self
        self._msg = msg
        super().__init__(msg)
    def GetType(self) -> str:
        return self._msg
    def GetAtomIdx(self) -> int:
        return getattr(self, "_atom_idx", None)
    def GetAtomIndices(self) -> list[int]:
        res = getattr(self, "_atom_indices", None)
        if res is not None:
          res = tuple(res)
        return res
)PY",
           scope);
  const auto exceptionBase = m.attr("_BaseSanitException");

  auto mse = nb::exception<MolSanitizeException>(m, "MolSanitizeException",
                                                 exceptionBase);

  auto ase =
      nb::exception<AtomSanitizeException>(m, "AtomSanitizeException", mse);
  // pattern for this from https://mlir.llvm.org/doxygen/IRCore_8cpp_source.html
  // it's gross that this needs to be done for every exception type,
  // but I did not find a better way to do it
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *payload) {
        try {
          if (p) std::rethrow_exception(p);
        } catch (AtomSanitizeException &e) {
          nb::object ty = nb::borrow(static_cast<PyObject *>(payload));
          nb::object obj = ty(e.what());
          obj.attr("_type") = nb::cast(e.getType());
          obj.attr("_atom_idx") = nb::cast(e.getAtomIdx());
          PyErr_SetObject(static_cast<PyObject *>(payload), obj.ptr());
        }
      },
      ase.ptr());
  auto ave =
      nb::exception<AtomValenceException>(m, "AtomValenceException", ase);
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *payload) {
        try {
          if (p) std::rethrow_exception(p);
        } catch (AtomValenceException &e) {
          nb::object ty = nb::borrow(static_cast<PyObject *>(payload));
          nb::object obj = ty(e.what());
          obj.attr("_type") = nb::cast(e.getType());
          obj.attr("_atom_idx") = nb::cast(e.getAtomIdx());
          PyErr_SetObject(static_cast<PyObject *>(payload), obj.ptr());
        }
      },
      ave.ptr());
  auto ake =
      nb::exception<AtomKekulizeException>(m, "AtomKekulizeException", ase);
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *payload) {
        try {
          if (p) std::rethrow_exception(p);
        } catch (AtomKekulizeException &e) {
          nb::object ty = nb::borrow(static_cast<PyObject *>(payload));
          nb::object obj = ty(e.what());
          obj.attr("_type") = nb::cast(e.getType());
          obj.attr("_atom_idx") = nb::cast(e.getAtomIdx());
          PyErr_SetObject(static_cast<PyObject *>(payload), obj.ptr());
        }
      },
      ake.ptr());
  auto ke = nb::exception<KekulizeException>(m, "KekulizeException", mse);
  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *payload) {
        try {
          if (p) std::rethrow_exception(p);
        } catch (KekulizeException &e) {
          nb::object ty = nb::borrow(static_cast<PyObject *>(payload));
          nb::object obj = ty(e.what());
          obj.attr("_type") = nb::cast(e.getType());
          obj.attr("_atom_indices") = nb::cast(e.getAtomIndices());
          PyErr_SetObject(static_cast<PyObject *>(payload), obj.ptr());
        }
      },
      ke.ptr());

  nb::class_<MolSanitizeException>(m, "_cppMolSanitizeException",
                                   "exception arising from sanitization")
      .def("Message", &MolSanitizeException::what)
      .def("GetType", &MolSanitizeException::getType)
      .def_prop_ro("cause",
                   [](const MolSanitizeException &self) { return self; });
  nb::class_<AtomSanitizeException, MolSanitizeException>(
      m, "_cppAtomSanitizeException", "exception arising from sanitization")
      .def("GetAtomIdx", &AtomSanitizeException::getAtomIdx);
  nb::class_<AtomValenceException, AtomSanitizeException>(
      m, "_cppAtomValenceException", "exception arising from sanitization");
  nb::class_<AtomKekulizeException, AtomSanitizeException>(
      m, "_cppAtomKekulizeException", "exception arising from sanitization");
  nb::class_<KekulizeException, MolSanitizeException>(
      m, "_cppKekulizeException", "exception arising from sanitization")
      .def("GetAtomIndices", &KekulizeException::getAtomIndices);

  //*********************************************
  //
  //  Classes
  //
  //*********************************************
  wrap_table(m);
  wrap_atom(m);
  wrap_bond(m);
  wrap_mol(m);
  wrap_conformer(m);
  wrap_stereogroup(m);
  wrap_EditableMol(m);
  wrap_ringinfo(m);
  wrap_monomerinfo(m);
  wrap_resmolsupplier(m);
  wrap_molbundle(m);
  wrap_sgroup(m);
  wrap_chirality(m);
}
