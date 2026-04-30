//
//  Copyright (C) 2026 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <DataStructs/BitOps.h>
#include <DataStructs/BitVects.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

using EBV = ExplicitBitVect;
using SBV = SparseBitVect;

namespace {
std::string bytesToString(const nb::bytes &b) {
  return std::string(static_cast<const char *>(b.data()),
                     static_cast<size_t>(b.size()));
}

template <typename T>
nb::bytes BVToBinaryText(const T &bv) {
  std::string res = BitVectToBinaryText(bv);
  return nb::bytes(res.c_str(), res.length());
}
}  // namespace

template <typename T>
double SimilarityWrapper(const T &bv1, const nb::bytes &pkl,
                         double (*metric)(const T &, const T &),
                         bool returnDistance) {
  T bv2(bytesToString(pkl));
  return SimilarityWrapper(bv1, bv2, metric, returnDistance);
}

template <typename T>
double SimilarityWrapper(const T &bv1, const nb::bytes &pkl, double a, double b,
                         double (*metric)(const T &, const T &, double, double),
                         bool returnDistance) {
  T bv2(bytesToString(pkl));
  return SimilarityWrapper(bv1, bv2, a, b, metric, returnDistance);
}

template <typename T>
nb::list NeighborWrapper(const nb::iterable &queries, const nb::iterable &bvs,
                         double (*metric)(const T &, const T &)) {
  nb::list res;
  std::vector<const T *> bvPtrs;
  for (auto item : bvs) {
    bvPtrs.push_back(&nb::cast<const T &>(item));
  }

  for (auto query : queries) {
    const auto &bv1 = nb::cast<const T &>(query);
    double closest = -1.0;
    unsigned int nbr = 0;
    for (size_t j = 0; j < bvPtrs.size(); ++j) {
      auto sim = metric(bv1, *bvPtrs[j]);
      if (sim > closest) {
        closest = sim;
        nbr = static_cast<unsigned int>(j);
      }
    }
    res.append(nb::make_tuple(nbr, closest));
  }
  return res;
}

template <typename T>
nb::list BulkWrapper(const T &bv1, const nb::iterable &bvs,
                     double (*metric)(const T &, const T &),
                     bool returnDistance) {
  nb::list res;
  for (auto item : bvs) {
    const auto &bv2 = nb::cast<const T &>(item);
    auto sim = metric(bv1, bv2);
    if (returnDistance) {
      sim = 1.0 - sim;
    }
    res.append(sim);
  }
  return res;
}

template <typename T>
nb::list BulkWrapper(const T &bv1, const nb::iterable &bvs, double a, double b,
                     double (*metric)(const T &, const T &, double, double),
                     bool returnDistance) {
  nb::list res;
  for (auto item : bvs) {
    const auto &bv2 = nb::cast<const T &>(item);
    auto sim = metric(bv1, bv2, a, b);
    if (returnDistance) {
      sim = 1.0 - sim;
    }
    res.append(sim);
  }
  return res;
}

#define METRIC_DEFS(_metricname_)                                              \
  template <typename T1, typename T2>                                          \
  double _metricname_##_w(const T1 &bv1, const T2 &bv2, bool returnDistance) { \
    return SimilarityWrapper(bv1, bv2,                                         \
                             (double (*)(const T1 &, const T1 &))_metricname_, \
                             returnDistance);                                  \
  }                                                                            \
  template <typename T>                                                        \
  nb::list Bulk##_metricname_(const T &bv1, const nb::iterable &bvs,           \
                              bool returnDistance) {                           \
    return BulkWrapper(bv1, bvs,                                               \
                       (double (*)(const T &, const T &))_metricname_,         \
                       returnDistance);                                        \
  }                                                                            \
  template <typename T>                                                        \
  nb::list _metricname_##Neighbors(const nb::iterable &queries,                \
                                   const nb::iterable &bvs) {                  \
    return NeighborWrapper<T>(queries, bvs,                                    \
                              (double (*)(const T &, const T &))_metricname_); \
  }

METRIC_DEFS(TanimotoSimilarity)
METRIC_DEFS(CosineSimilarity)
METRIC_DEFS(KulczynskiSimilarity)
METRIC_DEFS(DiceSimilarity)
METRIC_DEFS(SokalSimilarity)
METRIC_DEFS(McConnaugheySimilarity)
METRIC_DEFS(AsymmetricSimilarity)
METRIC_DEFS(BraunBlanquetSimilarity)
METRIC_DEFS(RusselSimilarity)
METRIC_DEFS(RogotGoldbergSimilarity)
METRIC_DEFS(OnBitSimilarity)
METRIC_DEFS(AllBitSimilarity)

template <typename T1, typename T2>
double TverskySimilarity_w(const T1 &bv1, const T2 &bv2, double a, double b,
                           bool returnDistance) {
  return SimilarityWrapper(
      bv1, bv2, a, b,
      (double (*)(const T1 &, const T1 &, double, double))TverskySimilarity,
      returnDistance);
}

template <typename T>
nb::list BulkTverskySimilarity(const T &bv1, const nb::iterable &bvs, double a,
                               double b, bool returnDistance) {
  return BulkWrapper(
      bv1, bvs, a, b,
      (double (*)(const T &, const T &, double, double))TverskySimilarity,
      returnDistance);
}

template <typename T>
bool AllProbeBitsMatchBytes(const T &probe, const nb::bytes &ref) {
  return AllProbeBitsMatch(probe, bytesToString(ref));
}

#define DBL_DEF(_funcname_, _bulkname_, _help_)                          \
  {                                                                      \
    m.def(#_funcname_, (double (*)(const SBV &, const SBV &))_funcname_, \
          "v1"_a, "v2"_a);                                               \
    m.def(#_funcname_, (double (*)(const EBV &, const EBV &))_funcname_, \
          "v1"_a, "v2"_a, _help_);                                       \
    m.def(#_bulkname_,                                                   \
          (nb::list (*)(const SBV &, const nb::iterable &,               \
                        bool))_bulkname_<SBV>,                           \
          "v1"_a, "v2"_a, "returnDistance"_a = false);                   \
    m.def(#_bulkname_,                                                   \
          (nb::list (*)(const EBV &, const nb::iterable &,               \
                        bool))_bulkname_<EBV>,                           \
          "v1"_a, "v2"_a, "returnDistance"_a = false, _help_);           \
  }

#define BIG_DEF(_funcname_, _name_w_, _bulkname_, _help_)                      \
  {                                                                            \
    m.def(#_funcname_, (double (*)(const SBV &, const SBV &, bool))_name_w_,   \
          "bv1"_a, "bv2"_a, "returnDistance"_a = false);                       \
    m.def(#_funcname_, (double (*)(const EBV &, const EBV &, bool))_name_w_,   \
          "bv1"_a, "bv2"_a, "returnDistance"_a = false, _help_);               \
    m.def(#_funcname_,                                                         \
          (double (*)(const SBV &, const nb::bytes &, bool))_name_w_, "bv1"_a, \
          "pkl"_a, "returnDistance"_a = false);                                \
    m.def(#_funcname_,                                                         \
          (double (*)(const EBV &, const nb::bytes &, bool))_name_w_, "bv1"_a, \
          "pkl"_a, "returnDistance"_a = false, _help_);                        \
    m.def(#_bulkname_,                                                         \
          (nb::list (*)(const SBV &, const nb::iterable &,                     \
                        bool))_bulkname_<SBV>,                                 \
          "bv1"_a, "bvList"_a, "returnDistance"_a = false);                    \
    m.def(#_bulkname_,                                                         \
          (nb::list (*)(const EBV &, const nb::iterable &,                     \
                        bool))_bulkname_<EBV>,                                 \
          "bv1"_a, "bvList"_a, "returnDistance"_a = false, _help_);            \
    m.def(#_funcname_ "Neighbors", _funcname_##Neighbors<ExplicitBitVect>,     \
          "bvqueries"_a, "bvList"_a, _help_);                                  \
    m.def(#_funcname_ "Neighbors_sparse",                                      \
          _funcname_##Neighbors<SparseBitVect>, "bvqueries"_a, "bvList"_a,     \
          _help_);                                                             \
  }

struct BitOps_wrapper {
  static void wrap(nb::module_ &m) {
    BIG_DEF(TanimotoSimilarity, TanimotoSimilarity_w, BulkTanimotoSimilarity,
            R"DOC(B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2)))DOC");
    BIG_DEF(CosineSimilarity, CosineSimilarity_w, BulkCosineSimilarity,
            R"DOC(B(bv1&bv2) / sqrt(B(bv1) * B(bv2)))DOC");
    BIG_DEF(KulczynskiSimilarity, KulczynskiSimilarity_w,
            BulkKulczynskiSimilarity,
            R"DOC(B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2)))DOC");
    BIG_DEF(DiceSimilarity, DiceSimilarity_w, BulkDiceSimilarity,
            R"DOC(2*B(bv1&bv2) / (B(bv1) + B(bv2)))DOC");
    BIG_DEF(SokalSimilarity, SokalSimilarity_w, BulkSokalSimilarity,
            R"DOC(B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2)))DOC");
    BIG_DEF(
        McConnaugheySimilarity, McConnaugheySimilarity_w,
        BulkMcConnaugheySimilarity,
        R"DOC((B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2)))DOC");
    BIG_DEF(AsymmetricSimilarity, AsymmetricSimilarity_w,
            BulkAsymmetricSimilarity,
            R"DOC(B(bv1&bv2) / min(B(bv1),B(bv2)))DOC");
    BIG_DEF(BraunBlanquetSimilarity, BraunBlanquetSimilarity_w,
            BulkBraunBlanquetSimilarity,
            R"DOC(B(bv1&bv2) / max(B(bv1),B(bv2)))DOC");
    BIG_DEF(RusselSimilarity, RusselSimilarity_w, BulkRusselSimilarity,
            R"DOC(B(bv1&bv2) / B(bv1))DOC");
    BIG_DEF(RogotGoldbergSimilarity, RogotGoldbergSimilarity_w,
            BulkRogotGoldbergSimilarity, R"DOC(B(bv1&bv2) / B(bv1))DOC");

    {
      const char *help =
          R"DOC(B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)))DOC";
      m.def("TverskySimilarity",
            (double (*)(const SBV &, const SBV &, double, double,
                        bool))TverskySimilarity_w,
            "bv1"_a, "bv2"_a, "a"_a, "b"_a, "returnDistance"_a = false);
      m.def("TverskySimilarity",
            (double (*)(const EBV &, const EBV &, double, double,
                        bool))TverskySimilarity_w,
            "bv1"_a, "bv2"_a, "a"_a, "b"_a, "returnDistance"_a = false, help);
      m.def("TverskySimilarity",
            (double (*)(const SBV &, const nb::bytes &, double, double,
                        bool))TverskySimilarity_w,
            "bv1"_a, "pkl"_a, "a"_a, "b"_a, "returnDistance"_a = false);
      m.def("TverskySimilarity",
            (double (*)(const EBV &, const nb::bytes &, double, double,
                        bool))TverskySimilarity_w,
            "bv1"_a, "pkl"_a, "a"_a, "b"_a, "returnDistance"_a = false, help);
      m.def("BulkTverskySimilarity",
            (nb::list (*)(const SBV &, const nb::iterable &, double, double,
                          bool))BulkTverskySimilarity<SBV>,
            "bv1"_a, "bvList"_a, "a"_a, "b"_a, "returnDistance"_a = false);
      m.def("BulkTverskySimilarity",
            (nb::list (*)(const EBV &, const nb::iterable &, double, double,
                          bool))BulkTverskySimilarity<EBV>,
            "bv1"_a, "bvList"_a, "a"_a, "b"_a, "returnDistance"_a = false,
            help);
    }

    DBL_DEF(OnBitSimilarity, BulkOnBitSimilarity,
            R"DOC(B(bv1&bv2) / B(bv1|bv2))DOC");
    DBL_DEF(AllBitSimilarity, BulkAllBitSimilarity,
            R"DOC((B(bv1) - B(bv1^bv2)) / B(bv1))DOC");

    m.def("OnBitProjSimilarity",
          (DoubleVect (*)(const SBV &, const SBV &))OnBitProjSimilarity,
          "bv1"_a, "bv2"_a);
    m.def(
        "OnBitProjSimilarity",
        (DoubleVect (*)(const EBV &, const EBV &))OnBitProjSimilarity, "bv1"_a,
        "bv2"_a,
        R"DOC(Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2)))DOC");
    m.def("OffBitProjSimilarity",
          (DoubleVect (*)(const SBV &, const SBV &))OffBitProjSimilarity,
          "bv1"_a, "bv2"_a);
    m.def("OffBitProjSimilarity",
          (DoubleVect (*)(const EBV &, const EBV &))OffBitProjSimilarity,
          "bv1"_a, "bv2"_a);

    m.def("NumBitsInCommon", (int (*)(const SBV &, const SBV &))NumBitsInCommon,
          "bv1"_a, "bv2"_a);
    m.def(
        "NumBitsInCommon", (int (*)(const EBV &, const EBV &))NumBitsInCommon,
        "bv1"_a, "bv2"_a,
        R"DOC(Returns the total number of bits in common between the two bit vectors)DOC");
    m.def("OnBitsInCommon",
          (IntVect (*)(const SBV &, const SBV &))OnBitsInCommon, "bv1"_a,
          "bv2"_a);
    m.def(
        "OnBitsInCommon", (IntVect (*)(const EBV &, const EBV &))OnBitsInCommon,
        "bv1"_a, "bv2"_a,
        R"DOC(Returns the number of on bits in common between the two bit vectors)DOC");
    m.def("OffBitsInCommon",
          (IntVect (*)(const SBV &, const SBV &))OffBitsInCommon, "bv1"_a,
          "bv2"_a);
    m.def(
        "OffBitsInCommon",
        (IntVect (*)(const EBV &, const EBV &))OffBitsInCommon, "bv1"_a,
        "bv2"_a,
        R"DOC(Returns the number of off bits in common between the two bit vectors)DOC");

    m.def("FoldFingerprint",
          (SBV * (*)(const SBV &, unsigned int)) FoldFingerprint, "bv"_a,
          "foldFactor"_a = 2, nb::rv_policy::take_ownership);
    m.def(
        "FoldFingerprint",
        (EBV * (*)(const EBV &, unsigned int)) FoldFingerprint, "bv"_a,
        "foldFactor"_a = 2, nb::rv_policy::take_ownership,
        R"DOC(Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.)DOC");

    m.def("AllProbeBitsMatch",
          (bool (*)(const SBV &, const SBV &))AllProbeBitsMatch, "probe"_a,
          "ref"_a);
    m.def("AllProbeBitsMatch",
          (bool (*)(const EBV &, const EBV &))AllProbeBitsMatch, "probe"_a,
          "ref"_a);
    m.def("AllProbeBitsMatch", AllProbeBitsMatchBytes<SBV>, "probe"_a, "ref"_a);
    m.def(
        "AllProbeBitsMatch", AllProbeBitsMatchBytes<EBV>, "probe"_a, "ref"_a,
        R"DOC(Returns True if all bits in the first argument match all bits in the
vector defined by the pickle in the second argument.
)DOC");

    m.def("BitVectToText", (std::string (*)(const SBV &))BitVectToText,
          "bv1"_a);
    m.def(
        "BitVectToText", (std::string (*)(const EBV &))BitVectToText, "bv1"_a,
        R"DOC(Returns a string of zeros and ones representing the bit vector.)DOC");
    m.def("BitVectToFPSText", (std::string (*)(const SBV &))BitVectToFPSText,
          "bv1"_a);
    m.def("BitVectToFPSText", (std::string (*)(const EBV &))BitVectToFPSText,
          "bv1"_a,
          R"DOC(Returns an FPS string representing the bit vector.)DOC");
    m.def("BitVectToBinaryText", (nb::bytes (*)(const SBV &))BVToBinaryText,
          "bv"_a);
    m.def(
        "BitVectToBinaryText", (nb::bytes (*)(const EBV &))BVToBinaryText,
        "bv"_a,
        R"DOC(Returns a binary string (byte array) representing the bit vector.)DOC");
  }
};

void wrap_BitOps(nb::module_ &m) { BitOps_wrapper::wrap(m); }
