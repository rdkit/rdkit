// Copyright Rational Discovery LLC 2005.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/shared_ptr.hpp>

#include <vector>

namespace python = boost::python;

class DemoKlass {
public:
  explicit DemoKlass(int v) : val_(v) {};
  int getVal() const { return val_; };
private:
  int val_;
};
typedef boost::shared_ptr<DemoKlass> DemoKlassSPtr;
typedef std::vector<DemoKlass*> DemoKlassPtrVect;
typedef std::vector<DemoKlassSPtr> DemoKlassSPtrVect;

DemoKlass *buildPtr(int v) {
  return new DemoKlass(v);
}
DemoKlassSPtr buildSPtr(int v) {
  return DemoKlassSPtr(new DemoKlass(v));
}

DemoKlassPtrVect buildPtrVector(unsigned int sz){
  DemoKlassPtrVect res;
  for(unsigned int i=0;i<sz;i++){
    res.push_back(new DemoKlass(i));
  }
  return res;
}

DemoKlassSPtrVect buildSPtrVector(unsigned int sz){
  DemoKlassSPtrVect res;
  for(unsigned int i=0;i<sz;i++){
    res.push_back(DemoKlassSPtr(new DemoKlass(i)));
  }
  return res;
}


class DemoContainer {
public:
  typedef DemoKlassSPtrVect::iterator iterator;
  typedef DemoKlassSPtrVect::const_iterator const_iterator;
  explicit DemoContainer(unsigned int sz) {
    vect_ = buildSPtrVector(sz);
  }
  iterator begin() {
    return vect_.begin();
  }
  iterator end() {
    return vect_.end();
  }
  const_iterator begin() const {
    return vect_.begin();
  }
  const_iterator end() const {
    return vect_.end();
  }

private:
  DemoKlassSPtrVect vect_;
};



BOOST_PYTHON_MODULE(SPtrTestModule)
{
  python::class_<DemoKlass,DemoKlassSPtr >("DemoKlass","demo class",python::init<int>())
    .def("GetVal",&DemoKlass::getVal)
    ;

  python::def("buildPtr",buildPtr,python::return_value_policy<python::manage_new_object>());
  python::def("buildSPtr",buildSPtr);


  python::class_<DemoKlassPtrVect>("DemoKlassPtrVec")
    .def(python::vector_indexing_suite<DemoKlassPtrVect>())
    ;
  python::def("buildPtrVector",buildPtrVector);

  python::class_<DemoKlassSPtrVect>("DemoKlassSPtrVec")
    .def(python::vector_indexing_suite<DemoKlassSPtrVect,true>())
    ;
  python::def("buildSPtrVector",buildSPtrVector);


  python::class_<DemoContainer>("DemoContainer","demo container",python::init<unsigned int>())
    .def("__iter__",python::iterator<DemoContainer>())
    ;


}
