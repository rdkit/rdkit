//
//  Copyright (C) 2026 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitVectUtils.h>
#include <DataStructs/BitOps.h>
#include <RDGeneral/Exceptions.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;
using namespace nb::literals;

ExplicitBitVect *createFromBitString(const std::string &bits) {
  auto *res = new ExplicitBitVect(bits.size());
  FromBitString(*res, bits);
  return res;
}

ExplicitBitVect *createFromFPSText(const std::string &fps) {
  if (fps.size() % 2) {
    throw ValueErrorException(
        "input string must have an even number of characters");
  }
  auto *res = new ExplicitBitVect(fps.size() * 4);
  UpdateBitVectFromFPSText(*res, fps);
  return res;
}

ExplicitBitVect *createFromBinaryText(const nb::bytes &fps) {
  auto *res = new ExplicitBitVect(fps.size() * 8);
  auto txt = std::string(static_cast<const char *>(fps.data()), fps.size());
  UpdateBitVectFromBinaryText(*res, txt);
  return res;
}
ExplicitBitVect *createFromBinaryText(const std::string &fps) {
  auto *res = new ExplicitBitVect(fps.size() * 8);
  UpdateBitVectFromBinaryText(*res, fps);
  return res;
}

struct Utils_wrapper {
  static void wrap(nb::module_ &m) {
    m.def(
        "ConvertToExplicit", convertToExplicit, "sbv"_a,
        nb::rv_policy::take_ownership,
        R"DOC(Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector)DOC");
    m.def(
        "CreateFromBitString", createFromBitString, "bits"_a,
        nb::rv_policy::take_ownership,
        R"DOC(Creates an ExplicitBitVect from a bit string (string of 0s and 1s).)DOC");
    m.def("CreateFromFPSText", createFromFPSText, "fps"_a,
          nb::rv_policy::take_ownership,
          R"DOC(Creates an ExplicitBitVect from an FPS string.)DOC");
    m.def(
        "CreateFromBinaryText",
        nb::overload_cast<const nb::bytes &>(createFromBinaryText), "fps"_a,
        nb::rv_policy::take_ownership,
        R"DOC(Creates an ExplicitBitVect from a binary string (byte array).)DOC");
    m.def("CreateFromBinaryText",
          nb::overload_cast<const std::string &>(createFromBinaryText), "fps"_a,
          nb::rv_policy::take_ownership,
          R"DOC(Creates an ExplicitBitVect from a string (byte array).)DOC");

    m.def(
        "InitFromDaylightString",
        (void (*)(SparseBitVect &, const std::string &))FromDaylightString,
        "sbv"_a, "s"_a,
        R"DOC(Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.

   **Arguments**
     - bv: either a _SparseBitVect_ or an _ExplicitBitVect_
     - txt: a string with the Daylight encoding (this is the text that
    the Daylight tools put in the FP field of a TDT)

)DOC");
    m.def(
        "InitFromDaylightString",
        (void (*)(ExplicitBitVect &, const std::string &))FromDaylightString,
        "sbv"_a, "s"_a,
        R"DOC(Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.

   **Arguments**
     - bv: either a _SparseBitVect_ or an _ExplicitBitVect_
     - txt: a string with the Daylight encoding (this is the text that
    the Daylight tools put in the FP field of a TDT)

)DOC");
  }
};

void wrap_Utils(nb::module_ &m) { Utils_wrapper::wrap(m); }
