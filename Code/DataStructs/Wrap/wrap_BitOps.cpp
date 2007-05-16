// $Id$
//
//  Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include <boost/python.hpp>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>


namespace python = boost::python;

// sometimes I <heart> C++ syntax:
//   the following piece of grossness is, apparently, how
//   you cast a pointer to a function that takes default
//   arguments.  If you leave out the "=2" bit, you don't
//   get the default arg values.
//SBV *(*ff1)(const SBV&,int=2) = FoldFingerprint;
// this is a fix for VC++, which doesn't seem to like the "=2" passed in above
SBV *ff1(const SBV &bv1, int factor=2) {
  return FoldFingerprint(bv1,factor);
}

//ff1 = &FoldFingerprint(const SBV&, int=2);
//EBV *(*ff2)(const EBV&,int) = FoldFingerprint;
// same fix as above for VC 7.1
EBV *ff2(const EBV &ev1, int factor=2) {
  return FoldFingerprint(ev1,factor);
}


template <typename T>
double SimilarityWrapper(const T &bv1,const std::string &pkl,
                         const double (*metric)(const T &,const T &)){
  double res=0.0;
  T bv2(pkl);
  if(bv1.GetNumBits()>bv2.GetNumBits()){
    T *bv1tmp = FoldFingerprint(bv1,bv1.GetNumBits()/bv2.GetNumBits());
    res = metric(*bv1tmp,bv2);
    delete bv1tmp;
  } else if(bv2.GetNumBits()>bv1.GetNumBits()){
    T *bv2tmp = FoldFingerprint(bv2,bv2.GetNumBits()/bv1.GetNumBits());
    res = metric(bv1,*bv2tmp);
    delete bv2tmp;
  } else {
    res = metric(bv1,bv2);
  }
  return res;
}

template <typename T>
python::list BulkWrapper(const T &bv1,python::list bvs,
                         const double (*metric)(const T &,const T &)){
  python::list res;
  unsigned int nbvs=python::extract<unsigned int>(bvs.attr("__len__")());
  for(unsigned int i=0;i<nbvs;++i){
    double simVal;
    T bv2=python::extract<T>(bvs[i])();
    if(bv1.GetNumBits()>bv2.GetNumBits()){
      T *bv1tmp = FoldFingerprint(bv1,bv1.GetNumBits()/bv2.GetNumBits());
      simVal = metric(*bv1tmp,bv2);
      delete bv1tmp;
    } else if(bv2.GetNumBits()>bv1.GetNumBits()){
      T *bv2tmp = FoldFingerprint(bv2,bv2.GetNumBits()/bv1.GetNumBits());
      simVal = metric(bv1,*bv2tmp);
      delete bv2tmp;
    } else {
      simVal = metric(bv1,bv2);
    }
    res.append(simVal);
  }
  return res;
}

template <typename T>
double TanimotoSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))TanimotoSimilarity);
}
template <typename T>
python::list BulkTanimotoSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))TanimotoSimilarity);
}

template <typename T>
double CosineSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))CosineSimilarity);
}
template <typename T>
python::list BulkCosineSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))CosineSimilarity);
}

template <typename T>
double KulczynskiSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))KulczynskiSimilarity);
}
template <typename T>
python::list BulkKulczynskiSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))KulczynskiSimilarity);
}


template <typename T>
double DiceSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))DiceSimilarity);
}
template <typename T>
python::list BulkDiceSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))DiceSimilarity);
}


template <typename T>
double SokalSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))SokalSimilarity);
}
template <typename T>
python::list BulkSokalSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))SokalSimilarity);
}


template <typename T>
double McConnaugheySimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))McConnaugheySimilarity);
}
template <typename T>
python::list BulkMcConnaugheySimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))McConnaugheySimilarity);
}


template <typename T>
double AsymmetricSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))AsymmetricSimilarity);
}
template <typename T>
python::list BulkAsymmetricSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))AsymmetricSimilarity);
}


template <typename T>
double BraunBlanquetSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))BraunBlanquetSimilarity);
}
template <typename T>
python::list BulkBraunBlanquetSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))BraunBlanquetSimilarity);
}


template <typename T>
double RusselSimilarity_w(const T &bv1,const std::string &pkl){
  return SimilarityWrapper(bv1,pkl,
                           (const double (*)(const T&,const T&))RusselSimilarity);
}
template <typename T>
python::list BulkRusselSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))RusselSimilarity);
}

template <typename T>
python::list BulkOnBitSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))OnBitSimilarity);
}
template <typename T>
python::list BulkAllBitSimilarity(const T &bv1,python::list bvs){
  return BulkWrapper(bv1,bvs,
                     (const double (*)(const T&,const T&))AllBitSimilarity);
}


BOOST_PYTHON_FUNCTION_OVERLOADS(sbv_fold_overloads, ff1, 1, 2);
BOOST_PYTHON_FUNCTION_OVERLOADS(ebv_fold_overloads, ff2, 1, 2);

#define DBL_DEF(_funcname_,_bulkname_,_help_) \
  python::def( # _funcname_,(const double (*)(const SBV &,const SBV &))_funcname_); \
  python::def( # _funcname_,(const double (*)(const EBV &,const EBV &))_funcname_,_help_);\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::list))_bulkname_);\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::list))_bulkname_,_help_);
  

#define BIG_DEF(_funcname_,_name_w_,_bulkname_,_help_) \
  python::def( # _funcname_,(const double (*)(const SBV &,const SBV &))_funcname_); \
  python::def( # _funcname_,(const double (*)(const EBV &,const EBV &))_funcname_,_help_);\
  python::def( # _funcname_,(double (*)(const SBV &,const std::string &))_name_w_);\
  python::def( # _funcname_,(double (*)(const EBV &,const std::string &))_name_w_,_help_);\
  python::def( # _bulkname_,(python::list (*)(const SBV &,python::list))_bulkname_);\
  python::def( # _bulkname_,(python::list (*)(const EBV &,python::list))_bulkname_,_help_);

struct BitOps_wrapper {
  static void wrap(){
    BIG_DEF(TanimotoSimilarity,TanimotoSimilarity_w,BulkTanimotoSimilarity,
            "B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))")
    BIG_DEF(CosineSimilarity,CosineSimilarity_w,BulkCosineSimilarity,
            "B(bv1&bv2) / sqrt(B(bv1) * B(bv2))")
    BIG_DEF(KulczynskiSimilarity,KulczynskiSimilarity_w,BulkKulczynskiSimilarity,
            "B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))")
    BIG_DEF(DiceSimilarity,DiceSimilarity_w,BulkDiceSimilarity,
            "2*B(bv1&bv2) / (B(bv1) + B(bv2))")
    BIG_DEF(SokalSimilarity,SokalSimilarity_w,BulkSokalSimilarity,
            "B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))")

    BIG_DEF(McConnaugheySimilarity,McConnaugheySimilarity_w,BulkMcConnaugheySimilarity,
            "(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))")
    BIG_DEF(AsymmetricSimilarity,AsymmetricSimilarity_w,BulkAsymmetricSimilarity,
            "B(bv1&bv2) / min(B(bv1),B(bv2))")
    BIG_DEF(BraunBlanquetSimilarity,BraunBlanquetSimilarity_w,BulkBraunBlanquetSimilarity,
            "B(bv1&bv2) / max(B(bv1),B(bv2))")
    BIG_DEF(RusselSimilarity,RusselSimilarity_w,BulkRusselSimilarity,
            "B(bv1&bv2) / B(bv1)")

    DBL_DEF(OnBitSimilarity,BulkOnBitSimilarity,
            "B(bv1&bv2) / B(bv1|bv2)")
    DBL_DEF(AllBitSimilarity,BulkAllBitSimilarity,
            "(B(bv1) - B(bv1^bv2)) / B(bv1)")
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
                (const int (*)(const SBV&,const SBV&))NumBitsInCommon);
    python::def("NumBitsInCommon",
                (const int (*)(const EBV&,const EBV&))NumBitsInCommon,
                "Returns the total number of bits in common between the two bit vectors\n"
                );
    python::def("OnBitsInCommon",
                (IntVect (*)(const SBV&,const SBV&))OnBitsInCommon);
    python::def("OnBitsInCommon",
                (IntVect (*)(const EBV&,const EBV&))OnBitsInCommon,
                "Returns the number of on bits in common between the two bit vectors\n"
                );
    python::def("OffBitsInCommon",
                (IntVect (*)(const SBV&,const SBV&))OffBitsInCommon);
    python::def("OffBitsInCommon",
                (IntVect (*)(const EBV&,const EBV&))OffBitsInCommon,
                "Returns the number of off bits in common between the two bit vectors\n"
                );

    python::def("FoldFingerprint",
                ff1,
                sbv_fold_overloads(python::args("bv","foldFactor")
                                   )[python::return_value_policy<python::manage_new_object>()]);
    python::def("FoldFingerprint",
                ff2,
                ebv_fold_overloads(python::args("bv","foldFactor"),
                                   "Returns a folded version of the bit vector\n"
                                   )[python::return_value_policy<python::manage_new_object>()]);

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

    python::def("BulkTanimotoSimilarity",(python::list (*)(const SBV &,python::list))BulkTanimotoSimilarity);
    python::def("BulkTanimotoSimilarity",(python::list (*)(const EBV &,python::list))BulkTanimotoSimilarity,
            "returns a list of similarity values");
  }
};

void wrap_BitOps() {
  BitOps_wrapper::wrap();
}
  
