//
//  Copyright (C) 2003-2021 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

namespace python = boost::python;

namespace {
template <typename T>
python::object BVToBinaryText(const T &bv) {
  std::string res = BitVectToBinaryText(bv);
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}
}  // namespace

template <typename T>
double SimilarityWrapper(const T &bv1, const std::string &pkl,
                         double (*metric)(const T &, const T &),
                         bool returnDistance) {
  T bv2(pkl);
  return SimilarityWrapper(bv1, bv2, metric, returnDistance);
}
template <typename T>
double SimilarityWrapper(const T &bv1, const std::string &pkl, double a,
                         double b,
                         double (*metric)(const T &, const T &, double, double),
                         bool returnDistance) {
  T bv2(pkl);
  return SimilarityWrapper(bv1, bv2, a, b, metric, returnDistance);
}

template <typename T>
python::list NeighborWrapper(python::object queries, python::object bvs,
                             double (*metric)(const T &, const T &)) {
  python::list res;
  unsigned int nbvs = python::extract<unsigned int>(bvs.attr("__len__")());
  unsigned int nqs = python::extract<unsigned int>(queries.attr("__len__")());
  for (unsigned int i = 0; i < nqs; ++i) {
    const T *bv1 = python::extract<const T *>(queries[i])();
    double closest = -1;
    unsigned nbr;
    for (unsigned int j = 0; j < nbvs; ++j) {
      const T *bv2 = python::extract<const T *>(bvs[j])();
      auto sim = metric(*bv1, *bv2);
      if (sim > closest) {
        closest = sim;
        nbr = j;
      }
    }
    res.append(python::make_tuple(nbr, closest));
  }
  return res;
}

template <typename T>
python::list BulkWrapper(const T *bv1, python::object bvs,
                         double (*metric)(const T &, const T &),
                         bool returnDistance) {
  python::list res;
  unsigned int nbvs = python::extract<unsigned int>(bvs.attr("__len__")());
  for (unsigned int i = 0; i < nbvs; ++i) {
    const T *bv2 = python::extract<const T *>(bvs[i])();
    auto sim = metric(*bv1, *bv2);
    if (returnDistance) {
      sim = 1.0 - sim;
    }
    res.append(sim);
  }
  return res;
}

template <typename T>
python::list BulkWrapper(const T *bv1, python::object bvs, double a, double b,
                         double (*metric)(const T &, const T &, double, double),
                         bool returnDistance) {
  python::list res;
  unsigned int nbvs = python::extract<unsigned int>(bvs.attr("__len__")());
  for (unsigned int i = 0; i < nbvs; ++i) {
    const T *bv2 = python::extract<T *>(bvs[i])();
    auto sim = metric(*bv1, *bv2, a, b);
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
  python::list Bulk##_metricname_(const T *bv1, python::object bvs,            \
                                  bool returnDistance) {                       \
    return BulkWrapper(bv1, bvs,                                               \
                       (double (*)(const T &, const T &))_metricname_,         \
                       returnDistance);                                        \
  }                                                                            \
  template <typename T>                                                        \
  python::list _metricname_##Neighbors(python::object queries,                 \
                                       python::object bvs) {                   \
    return NeighborWrapper(queries, bvs,                                       \
                           (double (*)(const T &, const T &))_metricname_);    \
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
python::list BulkTverskySimilarity(const T *bv1, python::object bvs, double a,
                                   double b, bool returnDistance) {
  return BulkWrapper(
      bv1, bvs, a, b,
      (double (*)(const T &, const T &, double, double))TverskySimilarity,
      returnDistance);
}

#define DBL_DEF(_funcname_, _bulkname_, _help_)                                \
  {                                                                            \
    python::def(#_funcname_, (double (*)(const SBV &, const SBV &))_funcname_, \
                (python::args("v1"), python::args("v2")));                     \
    python::def(#_funcname_, (double (*)(const EBV &, const EBV &))_funcname_, \
                (python::args("v1"), python::args("v2")), _help_);             \
    python::def(                                                               \
        #_bulkname_,                                                           \
        (python::list(*)(const EBV *, python::object, bool))_bulkname_,        \
        (python::args("v1"), python::args("v2"),                               \
         python::args("returnDistance") = 0));                                 \
    python::def(                                                               \
        #_bulkname_,                                                           \
        (python::list(*)(const EBV *, python::object, bool))_bulkname_,        \
        (python::args("v1"), python::args("v2"),                               \
         python::args("returnDistance") = 0),                                  \
        _help_);                                                               \
  }

#define BIG_DEF(_funcname_, _name_w_, _bulkname_, _help_)                     \
  {                                                                           \
    python::def(#_funcname_,                                                  \
                (double (*)(const SBV &, const SBV &, bool))_name_w_,         \
                (python::args("bv1"), python::args("bv2"),                    \
                 python::args("returnDistance") = 0));                        \
    python::def(#_funcname_,                                                  \
                (double (*)(const EBV &, const EBV &, bool))_name_w_,         \
                (python::args("bv1"), python::args("bv2"),                    \
                 python::args("returnDistance") = 0),                         \
                _help_);                                                      \
    python::def(#_funcname_,                                                  \
                (double (*)(const SBV &, const std::string &, bool))_name_w_, \
                (python::args("bv1"), python::args("pkl"),                    \
                 python::args("returnDistance") = 0));                        \
    python::def(#_funcname_,                                                  \
                (double (*)(const EBV &, const std::string &, bool))_name_w_, \
                (python::args("bv1"), python::args("pkl"),                    \
                 python::args("returnDistance") = 0),                         \
                _help_);                                                      \
    python::def(                                                              \
        #_bulkname_,                                                          \
        (python::list(*)(const SBV *, python::object, bool))_bulkname_,       \
        (python::args("bv1"), python::args("bvList"),                         \
         python::args("returnDistance") = 0));                                \
    python::def(                                                              \
        #_bulkname_,                                                          \
        (python::list(*)(const EBV *, python::object, bool))_bulkname_,       \
        (python::args("bv1"), python::args("bvList"),                         \
         python::args("returnDistance") = 0),                                 \
        _help_);                                                              \
    python::def(#_funcname_ "Neighbors",                                      \
                (python::list(*)(python::object, python::object))             \
                    _funcname_##Neighbors<ExplicitBitVect>,                   \
                (python::args("bvqueries"), python::args("bvList")), _help_); \
    python::def(#_funcname_ "Neighbors_sparse",                               \
                (python::list(*)(python::object, python::object))             \
                    _funcname_##Neighbors<SparseBitVect>,                     \
                (python::args("bvqueries"), python::args("bvList")), _help_); \
  }

struct BitOps_wrapper {
  static void wrap() {
    BIG_DEF(TanimotoSimilarity, TanimotoSimilarity_w, BulkTanimotoSimilarity,
            "B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))");
    BIG_DEF(CosineSimilarity, CosineSimilarity_w, BulkCosineSimilarity,
            "B(bv1&bv2) / sqrt(B(bv1) * B(bv2))");
    BIG_DEF(KulczynskiSimilarity, KulczynskiSimilarity_w,
            BulkKulczynskiSimilarity,
            "B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))");
    BIG_DEF(DiceSimilarity, DiceSimilarity_w, BulkDiceSimilarity,
            "2*B(bv1&bv2) / (B(bv1) + B(bv2))");
    BIG_DEF(SokalSimilarity, SokalSimilarity_w, BulkSokalSimilarity,
            "B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))");

    BIG_DEF(
        McConnaugheySimilarity, McConnaugheySimilarity_w,
        BulkMcConnaugheySimilarity,
        "(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))");
    BIG_DEF(AsymmetricSimilarity, AsymmetricSimilarity_w,
            BulkAsymmetricSimilarity, "B(bv1&bv2) / min(B(bv1),B(bv2))");
    BIG_DEF(BraunBlanquetSimilarity, BraunBlanquetSimilarity_w,
            BulkBraunBlanquetSimilarity, "B(bv1&bv2) / max(B(bv1),B(bv2))");
    BIG_DEF(RusselSimilarity, RusselSimilarity_w, BulkRusselSimilarity,
            "B(bv1&bv2) / B(bv1)");
    BIG_DEF(RogotGoldbergSimilarity, RogotGoldbergSimilarity_w,
            BulkRogotGoldbergSimilarity, "B(bv1&bv2) / B(bv1)");

    {
      std::string help = "B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)";
      python::def("TverskySimilarity",
                  (double (*)(const SBV &, const SBV &, double, double,
                              bool))TverskySimilarity_w,
                  (python::args("bv1"), python::args("bv2"), python::args("a"),
                   python::args("b"), python::args("returnDistance") = 0));
      python::def("TverskySimilarity",
                  (double (*)(const EBV &, const EBV &, double, double,
                              bool))TverskySimilarity_w,
                  (python::args("bv1"), python::args("bv2"), python::args("a"),
                   python::args("b"), python::args("returnDistance") = 0),
                  help.c_str());
      python::def("TverskySimilarity",
                  (double (*)(const SBV &, const std::string &, double, double,
                              bool))TverskySimilarity_w,
                  (python::args("bv1"), python::args("pkl"), python::args("a"),
                   python::args("b"), python::args("returnDistance") = 0));
      python::def("TverskySimilarity",
                  (double (*)(const EBV &, const std::string &, double, double,
                              bool))TverskySimilarity_w,
                  (python::args("bv1"), python::args("pkl"), python::args("a"),
                   python::args("b"), python::args("returnDistance") = 0),
                  help.c_str());
      python::def(
          "BulkTverskySimilarity",
          (python::list(*)(const SBV *, python::object, double, double,
                           bool))BulkTverskySimilarity,
          (python::args("bv1"), python::args("bvList"), python::args("a"),
           python::args("b"), python::args("returnDistance") = 0));
      python::def(
          "BulkTverskySimilarity",
          (python::list(*)(const EBV *, python::object, double, double,
                           bool))BulkTverskySimilarity,
          (python::args("bv1"), python::args("bvList"), python::args("a"),
           python::args("b"), python::args("returnDistance") = 0),
          help.c_str());
    }

    DBL_DEF(OnBitSimilarity, BulkOnBitSimilarity, "B(bv1&bv2) / B(bv1|bv2)");
    DBL_DEF(AllBitSimilarity, BulkAllBitSimilarity,
            "(B(bv1) - B(bv1^bv2)) / B(bv1)");

    python::def("OnBitProjSimilarity",
                (DoubleVect(*)(const SBV &, const SBV &))OnBitProjSimilarity);
    python::def(
        "OnBitProjSimilarity",
        (DoubleVect(*)(const EBV &, const EBV &))OnBitProjSimilarity,
        "Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))");
    python::def("OffBitProjSimilarity",
                (DoubleVect(*)(const SBV &, const SBV &))OffBitProjSimilarity);
    python::def("OffBitProjSimilarity",
                (DoubleVect(*)(const EBV &, const EBV &))OffBitProjSimilarity);

    python::def("NumBitsInCommon",
                (int (*)(const SBV &, const SBV &))NumBitsInCommon);
    python::def("NumBitsInCommon",
                (int (*)(const EBV &, const EBV &))NumBitsInCommon,
                "Returns the total number of bits in common between the two "
                "bit vectors");
    python::def("OnBitsInCommon",
                (IntVect(*)(const SBV &, const SBV &))OnBitsInCommon);
    python::def(
        "OnBitsInCommon", (IntVect(*)(const EBV &, const EBV &))OnBitsInCommon,
        "Returns the number of on bits in common between the two bit vectors");
    python::def("OffBitsInCommon",
                (IntVect(*)(const SBV &, const SBV &))OffBitsInCommon);
    python::def(
        "OffBitsInCommon",
        (IntVect(*)(const EBV &, const EBV &))OffBitsInCommon,
        "Returns the number of off bits in common between the two bit vectors");

    python::def("FoldFingerprint",
                (SBV * (*)(const SBV &, unsigned int)) FoldFingerprint,
                (python::arg("bv"), python::arg("foldFactor") = 2),
                python::return_value_policy<python::manage_new_object>());
    python::def("FoldFingerprint",
                (EBV * (*)(const EBV &, unsigned int)) FoldFingerprint,
                (python::arg("bv"), python::arg("foldFactor") = 2),
                python::return_value_policy<python::manage_new_object>(),
                "Folds the fingerprint by the provided amount. The default, "
                "foldFactor=2, returns a fingerprint that is half the size of "
                "the original.");

    python::def("AllProbeBitsMatch",
                (bool (*)(const SBV &, const SBV &))AllProbeBitsMatch);
    python::def("AllProbeBitsMatch",
                (bool (*)(const EBV &, const EBV &))AllProbeBitsMatch);
    python::def("AllProbeBitsMatch",
                (bool (*)(const SBV &, const std::string &))AllProbeBitsMatch);
    python::def(
        "AllProbeBitsMatch",
        (bool (*)(const EBV &, const std::string &))AllProbeBitsMatch,
        "Returns True if all bits in the first argument match all bits in the \n\
  vector defined by the pickle in the second argument.\n");

    python::def("BitVectToText", (std::string(*)(const SBV &))BitVectToText);
    python::def(
        "BitVectToText", (std::string(*)(const EBV &))BitVectToText,
        "Returns a string of zeros and ones representing the bit vector.");
    python::def("BitVectToFPSText",
                (std::string(*)(const SBV &))BitVectToFPSText);
    python::def("BitVectToFPSText",
                (std::string(*)(const EBV &))BitVectToFPSText,
                "Returns an FPS string representing the bit vector.");
    python::def("BitVectToBinaryText",
                (python::object(*)(const SBV &))BVToBinaryText);
    python::def(
        "BitVectToBinaryText", (python::object(*)(const EBV &))BVToBinaryText,
        "Returns a binary string (byte array) representing the bit vector.");
  }
};

void wrap_BitOps() { BitOps_wrapper::wrap(); }
