//
//  Copyright (C) 2026 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <DataStructs/BitVects.h>
#include <DataStructs/base64.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

using SBV = SparseBitVect;

void SetBitsFromList(SBV *bv, const nb::iterable &onBitList) {
  for (auto item : onBitList) {
    bv->setBit(nb::cast<unsigned int>(item));
  }
}

void UnSetBitsFromList(SBV *bv, const nb::iterable &offBitList) {
  for (auto item : offBitList) {
    bv->unsetBit(nb::cast<unsigned int>(item));
  }
}

int get_VectItem(const SBV &self, int which) {
  if (which < 0) {
    if (which + static_cast<int>(self.getNumBits()) < 0) {
      throw nb::index_error("bit index out of range");
    }
    which += self.getNumBits();
  }
  return self.getBit(static_cast<unsigned int>(which));
}

int set_VectItem(SBV &self, int which, const int val) {
  if (which < 0) {
    if (which + static_cast<int>(self.getNumBits()) < 0) {
      throw nb::index_error("bit index out of range");
    }
    which += self.getNumBits();
  }
  if (val) {
    return self.setBit(static_cast<unsigned int>(which));
  }
  return self.unsetBit(static_cast<unsigned int>(which));
}

IntVect GetOnBits(const SBV &self) {
  IntVect res;
  self.getOnBits(res);
  return res;
}

nb::bytes BVToBinary(const SBV &bv) {
  std::string res = bv.toString();
  return nb::bytes(res.c_str(), res.length());
}

void InitFromBase64(SBV &self, const std::string &inD) {
  self.initFromText(inD.c_str(), inD.length(), true);
}

std::string ToBase64(SBV &self) {
  std::string tmp = self.toString();
  const char *txt = Base64Encode(tmp.c_str(), tmp.length());
  std::string res(txt);
  delete[] txt;
  return res;
}

nb::list SparseToList(const SparseBitVect &sv) {
  nb::list l;
  if (sv.getNumBits()) {
    for (unsigned int i = 0; i < sv.getNumBits(); ++i) {
      l.append(0);
    }
    for (int i : *sv.getBitSet()) {
      l[static_cast<size_t>(i)] = 1;
    }
  }
  return l;
}

std::string sbvClassDoc =
    R"DOC(A class to store sparse bit vectors.

This class is most useful for situations where the size of the vector
is large and relatively few bits are set

For smaller or denser vectors, the _ExplicitBitVect_ class is much faster.

As you would expect, _SparseBitVects_ support a set of binary operations
so you can do things like:
     bv3 = bv1 & bv2  (bitwise and)
     bv3 = bv1 | bv2  (bitwise or)
     bv3 = bv1 ^ bv2  (bitwise xor)
     bv3 = ~bv1       (bitwise negation) NOTE: this operation is likely
                                                  to be VERY slow and inefficient.

Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
or by indexing (i.e. bv[i] = 1 or if bv[i]).

)DOC";
struct SBV_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<SparseBitVect>(m, "SparseBitVect", sbvClassDoc.c_str())
        .def(nb::init<unsigned int>(), "size"_a)
        .def(
            "__init__",
            [](SparseBitVect *t, nb::bytes b) {
              new (t)
                  SparseBitVect(std::string(static_cast<const char *>(b.data()),
                                            static_cast<size_t>(b.size())));
            },
            "pkl"_a)
        .def(
            "SetBit", (bool (SBV::*)(unsigned int))&SBV::setBit, "which"_a,
            R"DOC(Turns on a particular bit. Returns the original state of the bit.
)DOC")
        .def(
            "SetBitsFromList", SetBitsFromList, "onBitList"_a,
            R"DOC(Turns on a set of bits. The argument should be a tuple or list of bit ids.
)DOC")
        .def(
            "UnSetBit", (bool (SBV::*)(unsigned int))&SBV::unsetBit, "which"_a,
            R"DOC(Turns off a particular bit. Returns the original state of the bit.
)DOC")
        .def(
            "UnSetBitsFromList", UnSetBitsFromList, "offBitList"_a,
            R"DOC(Turns off a set of bits. The argument should be a tuple or list of bit ids.
)DOC")
        .def("GetBit", (bool (SBV::*)(unsigned int) const) & SBV::getBit,
             "which"_a, R"DOC(Returns the value of a bit.
)DOC")
        .def("GetNumBits", &SBV::getNumBits,
             R"DOC(Returns the number of bits in the vector (the vector's size).
)DOC")
        .def("__len__", &SBV::getNumBits)
        .def("GetNumOnBits", &SBV::getNumOnBits,
             R"DOC(Returns the number of on bits.
)DOC")
        .def("GetNumOffBits", &SBV::getNumOffBits,
             R"DOC(Returns the number of off bits.
)DOC")
        .def("__getitem__", get_VectItem, "which"_a)
        .def("__setitem__", set_VectItem, "which"_a, "val"_a)
        .def("GetOnBits", GetOnBits,
             R"DOC(Returns a tuple containing IDs of the on bits.
)DOC")
        .def("ToBinary", BVToBinary,
             R"DOC(Returns an internal binary representation of the vector.
)DOC")
        .def("FromBase64", InitFromBase64, "inD"_a,
             R"DOC(Initializes the vector from a base64 encoded binary string.
)DOC")
        .def(
            "ToBase64", ToBase64,
            R"DOC(Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
)DOC")
        .def("ToList", SparseToList,
             R"DOC(Return the BitVector as a python list.)DOC")
        .def(nb::self & nb::self)
        .def(nb::self | nb::self)
        .def(nb::self ^ nb::self)
        .def(~nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__getstate__",
             [](const SBV &sbv) {
               const auto pkl = sbv.toString();
               return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
             })
        .def("__setstate__",
             [](SBV &sbv, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&sbv) SBV(pkl);
             })
        .doc() = sbvClassDoc.c_str();
  }
};

void wrap_SBV(nb::module_ &m) { SBV_wrapper::wrap(m); }
