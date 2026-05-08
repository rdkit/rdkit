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

using EBV = ExplicitBitVect;

void SetBitsFromList(EBV *bv, const nb::iterable &onBitList) {
  for (auto item : onBitList) {
    bv->setBit(nb::cast<unsigned int>(item));
  }
}

void UnSetBitsFromList(EBV *bv, const nb::iterable &offBitList) {
  for (auto item : offBitList) {
    bv->unsetBit(nb::cast<unsigned int>(item));
  }
}

int get_VectItem(const EBV &self, int which) {
  if (which < 0) {
    if (which + static_cast<int>(self.getNumBits()) < 0) {
      throw nb::index_error("bit index out of range");
    }
    which += self.getNumBits();
  }
  return self.getBit(static_cast<unsigned int>(which));
}

int set_VectItem(EBV &self, int which, const int val) {
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

IntVect GetOnBits(const EBV &self) {
  IntVect res;
  self.getOnBits(res);
  return res;
}

nb::bytes BVToBinary(const EBV &bv) {
  std::string res = bv.toString();
  return nb::bytes(res.c_str(), res.length());
}

void InitFromBase64(EBV &self, const std::string &inD) {
  self.initFromText(inD.c_str(), inD.length(), true);
}

std::string ToBase64(EBV &self) {
  std::string tmp = self.toString();
  const char *txt = Base64Encode(tmp.c_str(), tmp.length());
  std::string res(txt);
  delete[] txt;
  return res;
}

nb::list ExplicitToList(const ExplicitBitVect &sv) {
  nb::list l;
  if (!sv.dp_bits) {
    return l;
  }

  auto count = sv.getNumBits();
  if (count) {
    for (unsigned int i = 0; i < count; ++i) {
      l.append(0);
    }
    auto pos = sv.dp_bits->find_first();
    while (pos != boost::dynamic_bitset<>::npos) {
      l[static_cast<size_t>(pos)] = 1;
      pos = sv.dp_bits->find_next(pos);
    }
  }
  return l;
}

std::string ebvClassDoc =
    R"DOC(A class to store explicit bit vectors.

This class is most useful for situations where the size of the vector
is relatively small (tens of thousands or smaller).

For larger vectors, use the _SparseBitVect_ class instead.

As you would expect, _ExplicitBitVects_ support a set of binary operations
so you can do things like:
     bv3 = bv1 & bv2  (bitwise and)
     bv3 = bv1 | bv2  (bitwise or)
     bv3 = bv1 ^ bv2  (bitwise xor)
     bv3 = ~bv1       (bitwise negation)

Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
or by indexing (i.e. bv[i] = 1 or if bv[i]).

)DOC";

struct EBV_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<ExplicitBitVect>(m, "ExplicitBitVect", ebvClassDoc.c_str())
        .def(nb::init<unsigned int>(), "size"_a)
        .def(
            "__init__",
            [](ExplicitBitVect *t, nb::bytes b) {
              new (t) ExplicitBitVect(
                  std::string(static_cast<const char *>(b.data()),
                              static_cast<size_t>(b.size())));
            },
            "pkl"_a)
        .def(nb::init<unsigned int, bool>(), "size"_a, "bitsSet"_a)
        .def(
            "SetBit", (bool (EBV::*)(unsigned int))&EBV::setBit, "which"_a,
            R"DOC(Turns on a particular bit. Returns the original state of the bit.
)DOC")
        .def(
            "SetBitsFromList", SetBitsFromList, "onBitList"_a,
            R"DOC(Turns on a set of bits. The argument should be a tuple or list of bit ids.
)DOC")
        .def(
            "UnSetBit", (bool (EBV::*)(unsigned int))&EBV::unsetBit, "which"_a,
            R"DOC(Turns off a particular bit. Returns the original state of the bit.
)DOC")
        .def(
            "UnSetBitsFromList", UnSetBitsFromList, "offBitList"_a,
            R"DOC(Turns off a set of bits. The argument should be a tuple or list of bit ids.
)DOC")
        .def("GetBit", (bool (EBV::*)(unsigned int) const) & EBV::getBit,
             "which"_a, R"DOC(Returns the value of a bit.
)DOC")
        .def("GetNumBits", &EBV::getNumBits,
             R"DOC(Returns the number of bits in the vector (the vector's size).
)DOC")
        .def("__len__", &EBV::getNumBits)
        .def("GetNumOnBits", &EBV::getNumOnBits,
             R"DOC(Returns the number of on bits.
)DOC")
        .def("GetNumOffBits", &EBV::getNumOffBits,
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
        .def(
            "ToList", ExplicitToList,
            R"DOC(Return the Bitvector as a python list (faster than list(vect)))DOC")
        .def(nb::self & nb::self)
        .def(nb::self | nb::self)
        .def(nb::self ^ nb::self)
        .def(nb::self + nb::self)
        .def(~nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self += nb::self)
        .def("__getstate__",
             [](const EBV &ebv) {
               const auto pkl = ebv.toString();
               return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
             })
        .def("__setstate__",
             [](EBV &ebv, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&ebv) EBV(pkl);
             })
        .doc() = ebvClassDoc.c_str();
  }
};

void wrap_EBV(nb::module_ &m) { EBV_wrapper::wrap(m); }
