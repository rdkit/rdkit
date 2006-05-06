// $Id: Queries.cpp 4978 2006-02-18 00:59:33Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit{


  template <class t1,class t2>
struct queries_wrapper {
  static void wrap(const char *type1Name,const char *type2Name){
    std::string className=type1Name;
    className += type2Name;
    python::class_<Queries::Query<t1,t2> >((className+"Query").c_str())
      .def("setNegation",&Queries::Query<t1,t2>::setNegation)
      ;

    python::class_< Queries::AndQuery<t1,t2> >((className+"AndQuery").c_str())
      ;
  };
};
}// end of namespace


void wrap_queries() {
  RDKit::queries_wrapper<int,int>::wrap("Int","Int");
}


