// $Id$
//
//  Copyright (C) 2007-2008 Greg Landrum
//
//  @@ All Rights Reserved @@
//
#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <DataStructs/SparseIntVect.h>
#include <boost/cstdint.hpp>

using namespace RDKit;

template <typename IndexType>
struct siv_pickle_suite : python::pickle_suite
{
  static python::tuple
  getinitargs(const SparseIntVect<IndexType>& self)
  {
    return python::make_tuple(self.toString());
  };
};

namespace {
  template <typename IndexType>
  void pyUpdateFromSequence(SparseIntVect<IndexType> &vect,
			    python::object &seq){
    PySequenceHolder<IndexType> seqL(seq);
    for(unsigned int i=0;i<seqL.size();++i){
      IndexType idx=seqL[i];
      vect.setVal(idx,vect[idx]+1);
    }
  }
  template <typename IndexType>
  python::dict pyGetNonzeroElements(SparseIntVect<IndexType> &vect){
    python::dict res;
    typename SparseIntVect<IndexType>::StorageType::const_iterator iter=vect.getNonzeroElements().begin();
    while(iter!=vect.getNonzeroElements().end()){
      res[iter->first]=iter->second;
      ++iter;
    }
    return res;
  }


  template <typename T>
  python::list BulkWrapper(const T &siv1,python::list sivs,bool returnDistance){
    python::list res;
    unsigned int nsivs=python::extract<unsigned int>(sivs.attr("__len__")());
    for(unsigned int i=0;i<nsivs;++i){
      double simVal;
      const T &siv2=python::extract<T>(sivs[i])();
      simVal = DiceSimilarity(siv1,siv2,returnDistance);
      res.append(simVal);
    }
    return res;
  }
}

std::string sparseIntVectDoc="A container class for storing integer\n\
values within a particular range.\n\
\n\
The length of the vector is set at construction time.\n\
\n\
As you would expect, _SparseIntVects_ support a set of binary operations\n\
so you can do things like:\n\
  Arithmetic:\n\
  siv1 += siv2\n\
  siv3 = siv1 + siv2\n\
  siv1 -= siv3\n\
  siv3 = siv1 - siv2\n\
  \"Fuzzy\" binary operations:\n\
  siv3 = siv1 & siv2  the result contains the smallest value in each entry\n\
  siv3 = siv1 | siv2  the result contains the largest value in each entry\n\
\n\
Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])\n\
\n";

struct sparseIntVec_wrapper {
  template <typename IndexType>
  static void wrapOne(const char *className){
    python::class_<SparseIntVect<IndexType> >(className, 
					      sparseIntVectDoc.c_str(),
					      python::init<IndexType>("Constructor"))
      .def(python::init<std::string>())
      // Note: we cannot support __len__ because, at least at the moment
      // (BPL v1.34.1), it must return an int.
      .def("__setitem__", &SparseIntVect<IndexType>::setVal,
	   "Set the value at a specified location")
      .def("__getitem__", &SparseIntVect<IndexType>::getVal,
	   "Get the value at a specified location")
      .def(python::self & python::self)
      .def(python::self | python::self)
      .def(python::self - python::self)
      .def(python::self -= python::self)
      .def(python::self + python::self)
      .def(python::self += python::self)
      .def(python::self == python::self)
      .def(python::self != python::self)
      .def("GetTotalVal", &SparseIntVect<IndexType>::getTotalVal,
	   (python::args("useAbs")=false),
	   "Get the sum of the values in the vector, basically L1 norm")
      .def("GetLength", &SparseIntVect<IndexType>::getLength,
	   (python::args("useAbs")=false),
	   "Returns the length of the vector")
      .def("ToBinary", &SparseIntVect<IndexType>::toString,
	   "returns a binary (pickle) representation of the vector")
      .def("UpdateFromSequence",
	   &pyUpdateFromSequence<IndexType>,
	   "update the vector based on the values in the list or tuple")
      .def("GetNonzeroElements",
	   &pyGetNonzeroElements<IndexType>,
	   "returns a dictionary of the nonzero elements")
      .def_pickle(siv_pickle_suite<IndexType>())
      ;

    python::def("DiceSimilarity",&DiceSimilarity<IndexType>,
		(python::args("siv1"),python::args("siv2"),
		 python::args("returnDistance")=false,
		 python::args("bounds")=0.0),
                "return the Dice similarity between two vectors");
    python::def("BulkDiceSimilarity",&BulkWrapper<SparseIntVect<IndexType> >,
		(python::args("v1"),python::args("v2"),
                 python::args("returnDistance")=false),
                "return the Dice similarities between one vector and a sequence of others");
  }

  static void wrap() {
    wrapOne<boost::int32_t>("IntSparseIntVect");
    wrapOne<boost::int64_t>("LongSparseIntVect");
  }
};
  
void wrap_sparseIntVect() {
  sparseIntVec_wrapper::wrap();
}

