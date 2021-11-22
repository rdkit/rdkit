// $Id$
//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//
// Generic Wrapper utility functionality
//
#include "Wrap.h"
#include "pyint_api.h"
#include <RDBoost/PySequenceHolder.h>
#include <sstream>
#include <iostream>

// A helper function for dealing with errors. Throw a Python IndexError
void throw_index_error(int key) {
  PyErr_SetObject(PyExc_IndexError, PyInt_FromLong(key));
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python ValueError
void throw_value_error(const std::string err) {
  PyErr_SetString(PyExc_ValueError, err.c_str());
  python::throw_error_already_set();
}

// A helper function for dealing with errors. Throw a Python KeyError
void throw_key_error(const std::string key) {
  PyErr_SetString(PyExc_KeyError, key.c_str());
  python::throw_error_already_set();
}

void translate_index_error(IndexErrorException const &e) {
  throw_index_error(e.index());
}

void translate_value_error(ValueErrorException const &e) {
  throw_value_error(e.what());
}

void translate_key_error(KeyErrorException const &e) {
  throw_key_error(e.key());
}

#ifdef INVARIANT_EXCEPTION_METHOD
// A helper function for dealing with errors. Throw a Python RuntimeError
void throw_runtime_error(const std::string err) {
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
  python::throw_error_already_set();
}

void translate_invariant_error(Invar::Invariant const &e) {
  throw_runtime_error(e.toUserString());
}

boost::dynamic_bitset<> pythonObjectToDynBitset(
    const python::object &obj, boost::dynamic_bitset<>::size_type maxV) {
  boost::dynamic_bitset<> res(maxV);
  if (obj) {
    python::stl_input_iterator<boost::dynamic_bitset<>::size_type> beg(obj),
        end;
    while (beg != end) {
      auto v = *beg;
      if (v >= maxV) {
        throw_value_error("list element larger than allowed value");
      }
      res.set(v);
      ++beg;
    }
  }
  return res;
}

std::vector<std::pair<int, int>> *translateAtomMap(
    const python::object &atomMap) {
  PySequenceHolder<python::object> pyAtomMap(atomMap);
  std::vector<std::pair<int, int>> *res;
  res = nullptr;
  unsigned int i;
  unsigned int n = pyAtomMap.size();
  if (n > 0) {
    res = new std::vector<std::pair<int, int>>;
    for (i = 0; i < n; ++i) {
      PySequenceHolder<int> item(pyAtomMap[i]);
      if (item.size() != 2) {
        delete res;
        res = nullptr;
        throw_value_error("Incorrect format for an atomMap");
      }
      res->push_back(std::pair<int, int>(item[0], item[1]));
    }
  }
  return res;
}

std::vector<std::vector<std::pair<int, int>>> translateAtomMapSeq(
    const python::object &atomMapSeq) {
  std::vector<std::vector<std::pair<int, int>>> aMapVec;
  PySequenceHolder<python::object> pyAtomMapSeq(atomMapSeq);
  for (size_t i = 0; i < pyAtomMapSeq.size(); ++i) {
    std::vector<std::pair<int, int>> *res = translateAtomMap(pyAtomMapSeq[i]);
    aMapVec.push_back(*res);
    delete res;
  }
  return aMapVec;
}

RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
  PySequenceHolder<double> doubles(doubleSeq);
  unsigned int nDoubles = doubles.size();
  RDNumeric::DoubleVector *doubleVec;
  doubleVec = nullptr;
  unsigned int i;
  if (nDoubles > 0) {
    doubleVec = new RDNumeric::DoubleVector(nDoubles);
    for (i = 0; i < nDoubles; ++i) {
      doubleVec->setVal(i, doubles[i]);
    }
  }
  return doubleVec;
}

std::vector<unsigned int> *translateIntSeq(const python::object &intSeq) {
  PySequenceHolder<unsigned int> ints(intSeq);
  std::vector<unsigned int> *intVec = nullptr;
  if (ints.size() > 0) {
    intVec = new std::vector<unsigned int>;
    for (unsigned int i = 0; i < ints.size(); ++i) {
      intVec->push_back(ints[i]);
    }
  }
  return intVec;
}

#endif
