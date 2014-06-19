// $Id$
//
//  Copyright (C) 2003-2012 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>


namespace python = boost::python;

SBV *ff1(const SBV &bv1, int factor=2) {
  return FoldFingerprint(bv1,factor);
}
EBV *ff2(const EBV &ev1, int factor=2) {
  return FoldFingerprint(ev1,factor);
}


namespace {
  template <typename T>
  python::object BVToBinaryText(const T &bv){
    std::string res=BitVectToBinaryText(bv);
    python::object retval = python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length())));
    return retval;
  }
}

template <typename T>
double SimilarityWrapper(const T &bv1,const std::string &pkl,
                         double (*metric)(const T &,const T &),bool returnDistance){
  T bv2(pkl);
  return SimilarityWrapper(bv1,bv2,metric,returnDistance);
}
template <typename T>
double SimilarityWrapper(const T &bv1,const std::string &pkl,double a,double b,
                         double (*metric)(const T &,const T &,double,double),bool returnDistance){
  T bv2(pkl);
  return SimilarityWrapper(bv1,bv2,a,b,metric,returnDistance);
}


template <typename T>
python::list BulkWrapper(const T &bv1,python::object bvs,
                         double (*metric)(const T &,const T &),
                         bool returnDistance){
  python::list res;
  unsigned int nbvs=python::extract<unsigned int>(bvs.attr("__len__")());
  for(unsigned int i=0;i<nbvs;++i){
    const T &bv2=python::extract<T>(bvs[i])();
    res.append(SimilarityWrapper(bv1,bv2,metric,returnDistance));
  }
  return res;
}

template <typename T>
python::list BulkWrapper(const T &bv1,python::object bvs,double a,double b,
                         double (*metric)(const T &,const T &,double,double),
                         bool returnDistance){
  python::list res;
  unsigned int nbvs=python::extract<unsigned int>(bvs.attr("__len__")());
  for(unsigned int i=0;i<nbvs;++i){
    const T &bv2=python::extract<T>(bvs[i])();
    res.append(SimilarityWrapper(bv1,bv2,a,b,metric,returnDistance));
  }
  return res;
}

template <typename T1, typename T2>
double TanimotoSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))TanimotoSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkTanimotoSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))TanimotoSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double TverskySimilarity_w(const T1 &bv1,const T2 &bv2,double a,double b,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,a,b,
                           (double (*)(const T1&,const T1&,double,double))TverskySimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkTverskySimilarity(const T &bv1,python::object bvs,double a,double b,
                                   bool returnDistance){
  return BulkWrapper(bv1,bvs,a,b,
                     (double (*)(const T&,const T&,double,double))TverskySimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double CosineSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))CosineSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkCosineSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))CosineSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double KulczynskiSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))KulczynskiSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkKulczynskiSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))KulczynskiSimilarity,
                     returnDistance);
}


template <typename T1, typename T2>
double DiceSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))DiceSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkDiceSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))DiceSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double SokalSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))SokalSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkSokalSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))SokalSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double McConnaugheySimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))McConnaugheySimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkMcConnaugheySimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))McConnaugheySimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double AsymmetricSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))AsymmetricSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkAsymmetricSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))AsymmetricSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double BraunBlanquetSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))BraunBlanquetSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkBraunBlanquetSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))BraunBlanquetSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double RusselSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))RusselSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkRusselSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))RusselSimilarity,
                     returnDistance);
}

template <typename T1, typename T2>
double RogotGoldbergSimilarity_w(const T1 &bv1,const T2 &bv2,bool returnDistance){
  return SimilarityWrapper(bv1,bv2,
                           (double (*)(const T1&,const T1&))RogotGoldbergSimilarity,
                           returnDistance);
}
template <typename T>
python::list BulkRogotGoldbergSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))RogotGoldbergSimilarity,
                     returnDistance);
}

template <typename T>
python::list BulkOnBitSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))OnBitSimilarity,
                     returnDistance);
}
template <typename T>
python::list BulkAllBitSimilarity(const T &bv1,python::object bvs,bool returnDistance){
  return BulkWrapper(bv1,bvs,
                     (double (*)(const T&,const T&))AllBitSimilarity,
                     returnDistance);
}


#define DBL_DEF(_funcname_,_bulkname_,_help_) { \
  python::def( # _funcname_,(double (*)(const SBV &,const SBV &))_funcname_,\
               (python::args("v1"),python::args("v2"))); \
  python::def( # _funcname_,(double (*)(const EBV &,const EBV &))_funcname_,\
               (python::args("v1"),python::args("v2")),_help_);\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::object,bool))_bulkname_,\
               (python::args("v1"),python::args("v2"),python::args("returnDistance")=0));\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::object,bool))_bulkname_,\
               (python::args("v1"),python::args("v2"),python::args("returnDistance")=0),_help_);}
  

#define BIG_DEF(_funcname_,_name_w_,_bulkname_,_help_) { \
  python::def( # _funcname_,(double (*)(const SBV &,const SBV &,bool))_name_w_,\
               (python::args("bv1"),python::args("bv2"),python::args("returnDistance")=0)); \
  python::def( # _funcname_,(double (*)(const EBV &,const EBV &,bool))_name_w_,\
               (python::args("bv1"),python::args("bv2"),python::args("returnDistance")=0),_help_);\
  python::def( # _funcname_,(double (*)(const SBV &,const std::string &,bool))_name_w_,\
               (python::args("bv1"),python::args("pkl"),python::args("returnDistance")=0));\
  python::def( # _funcname_,(double (*)(const EBV &,const std::string &,bool))_name_w_,\
               (python::args("bv1"),python::args("pkl"),python::args("returnDistance")=0),_help_);\
  python::def( # _bulkname_,(python::list (*)(const SBV &,python::object,bool))_bulkname_,\
               (python::args("bv1"),python::args("bvList"),python::args("returnDistance")=0));\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::object,bool))_bulkname_,\
               (python::args("bv1"),python::args("bvList"),python::args("returnDistance")=0),_help_);}

struct BitOps_wrapper {
  static void wrap(){
    BIG_DEF(TanimotoSimilarity,TanimotoSimilarity_w,BulkTanimotoSimilarity,
            "B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))");
    BIG_DEF(CosineSimilarity,CosineSimilarity_w,BulkCosineSimilarity,
            "B(bv1&bv2) / sqrt(B(bv1) * B(bv2))");
    BIG_DEF(KulczynskiSimilarity,KulczynskiSimilarity_w,BulkKulczynskiSimilarity,
            "B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))");
    BIG_DEF(DiceSimilarity,DiceSimilarity_w,BulkDiceSimilarity,
            "2*B(bv1&bv2) / (B(bv1) + B(bv2))");
    BIG_DEF(SokalSimilarity,SokalSimilarity_w,BulkSokalSimilarity,
            "B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))");

    BIG_DEF(McConnaugheySimilarity,McConnaugheySimilarity_w,BulkMcConnaugheySimilarity,
            "(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))");
    BIG_DEF(AsymmetricSimilarity,AsymmetricSimilarity_w,BulkAsymmetricSimilarity,
            "B(bv1&bv2) / min(B(bv1),B(bv2))");
    BIG_DEF(BraunBlanquetSimilarity,BraunBlanquetSimilarity_w,BulkBraunBlanquetSimilarity,
            "B(bv1&bv2) / max(B(bv1),B(bv2))");
    BIG_DEF(RusselSimilarity,RusselSimilarity_w,BulkRusselSimilarity,
            "B(bv1&bv2) / B(bv1)");
    BIG_DEF(RogotGoldbergSimilarity,RogotGoldbergSimilarity_w,BulkRogotGoldbergSimilarity,
                "B(bv1&bv2) / B(bv1)");

    {
      std::string help="B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)";
      python::def( "TverskySimilarity",
                   (double (*)(const SBV &,const SBV &,double,double,bool))TverskySimilarity_w,
                   (python::args("bv1"),python::args("bv2"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0));
      python::def( "TverskySimilarity",
                   (double (*)(const EBV &,const EBV &,double,double,bool))TverskySimilarity_w,
                   (python::args("bv1"),python::args("bv2"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0),help.c_str());
      python::def( "TverskySimilarity",
                   (double (*)(const SBV &,const std::string &,double,double,bool))TverskySimilarity_w,
                   (python::args("bv1"),python::args("pkl"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0));
      python::def( "TverskySimilarity",
                   (double (*)(const EBV &,const std::string &,double,double,bool))TverskySimilarity_w,
                   (python::args("bv1"),python::args("pkl"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0),help.c_str());
      python::def( "BulkTverskySimilarity",
                   (python::list (*)(const SBV &,python::object,double,double,bool))BulkTverskySimilarity,
                   (python::args("bv1"),python::args("bvList"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0));
      python::def( "BulkTverskySimilarity",
                   (python::list (*)(const EBV &,python::object,double,double,bool))BulkTverskySimilarity,
                   (python::args("bv1"),python::args("bvList"),python::args("a"),
                    python::args("b"),python::args("returnDistance")=0),help.c_str());
    }

    DBL_DEF(OnBitSimilarity,BulkOnBitSimilarity,
            "B(bv1&bv2) / B(bv1|bv2)");
    DBL_DEF(AllBitSimilarity,BulkAllBitSimilarity,
            "(B(bv1) - B(bv1^bv2)) / B(bv1)");


    python::def("OnBitProjSimilarity",
                (DoubleVect (*)(const SBV&,const SBV&))OnBitProjSimilarity);
    python::def("OnBitProjSimilarity",
                (DoubleVect (*)(const EBV&,const EBV&))OnBitProjSimilarity,
                "Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))");
    python::def("OffBitProjSimilarity",
                (DoubleVect (*)(const SBV&,const SBV&))OffBitProjSimilarity);
    python::def("OffBitProjSimilarity",
                (DoubleVect (*)(const EBV&,const EBV&))OffBitProjSimilarity);

    python::def("NumBitsInCommon",
                (int (*)(const SBV&,const SBV&))NumBitsInCommon);
    python::def("NumBitsInCommon",
                (int (*)(const EBV&,const EBV&))NumBitsInCommon,
                "Returns the total number of bits in common between the two bit vectors"
                );
    python::def("OnBitsInCommon",
                (IntVect (*)(const SBV&,const SBV&))OnBitsInCommon);
    python::def("OnBitsInCommon",
                (IntVect (*)(const EBV&,const EBV&))OnBitsInCommon,
                "Returns the number of on bits in common between the two bit vectors"
                );
    python::def("OffBitsInCommon",
                (IntVect (*)(const SBV&,const SBV&))OffBitsInCommon);
    python::def("OffBitsInCommon",
                (IntVect (*)(const EBV&,const EBV&))OffBitsInCommon,
                "Returns the number of off bits in common between the two bit vectors"
                );

    python::def("FoldFingerprint",
                (SBV *(*)(const SBV &,unsigned int))FoldFingerprint,
                (python::arg("bv"),python::arg("foldFactor")=2),
                python::return_value_policy<python::manage_new_object>());
    python::def("FoldFingerprint",
                (EBV *(*)(const EBV &,unsigned int))FoldFingerprint,
                (python::arg("bv"),python::arg("foldFactor")=2),
                python::return_value_policy<python::manage_new_object>(),
                "Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.");


    python::def("AllProbeBitsMatch",
                (bool (*)(const SBV &,const SBV &))AllProbeBitsMatch);
    python::def("AllProbeBitsMatch",
                (bool (*)(const EBV &,const EBV &))AllProbeBitsMatch);
    python::def("AllProbeBitsMatch",
                (bool (*)(const SBV &,const std::string &))AllProbeBitsMatch);
    python::def("AllProbeBitsMatch",(bool (*)(const EBV &,const std::string &))AllProbeBitsMatch,
                "Returns True if all bits in the first argument match all bits in the \n\
  vector defined by the pickle in the second argument.\n");

    python::def("BitVectToText",
                (std::string (*)(const SBV&))BitVectToText);
    python::def("BitVectToText",
                (std::string (*)(const EBV&))BitVectToText,
                "Returns a string of zeros and ones representing the bit vector."
                );
    python::def("BitVectToFPSText",
                (std::string (*)(const SBV&))BitVectToFPSText);
    python::def("BitVectToFPSText",
                (std::string (*)(const EBV&))BitVectToFPSText,
                "Returns an FPS string representing the bit vector."
                );
    python::def("BitVectToBinaryText",
                (python::object (*)(const SBV&))BVToBinaryText);
    python::def("BitVectToBinaryText",
                (python::object (*)(const EBV&))BVToBinaryText,
                "Returns a binary string (byte array) representing the bit vector."
                );
  }
};

void wrap_BitOps() {
  BitOps_wrapper::wrap();
}
  
