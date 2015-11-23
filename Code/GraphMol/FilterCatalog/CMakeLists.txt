INCLUDE (CheckCXXSourceCompiles)

set(T_LIBS ${Boost_LIBRARIES})
find_package(Boost 1.39.0 COMPONENTS serialization)
if (Boost_SERIALIZATION_LIBRARY)
  set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
  message("== Using boost serialization library ${Boost_SERIALIZATION_LIBRARY}")
  ADD_DEFINITIONS("-DRDK_USE_BOOST_SERIALIZATION")
else()
  message("== Making FilterCatalog without boost Serialization support")
endif()

# check to see if shared<const T> is supported by boost
IF(WIN32)
   add_definitions("-DBOOST_PYTHON_SUPPORT_SHARED_CONST")
else()
  set(CMAKE_REQUIRED_INCLUDES "${PYTHON_INCLUDE_DIR}")
  set(CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES};${Boost_INCLUDE_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES "${Boost_LIBRARIES};${PYTHON_LIBRARY}")
  CHECK_CXX_SOURCE_COMPILES("#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/register_ptr_to_python.hpp>
namespace python = boost::python;
namespace {
class Foo {
public:
  int a;
};
struct wrapper {
  static void wrap() {
      python::register_ptr_to_python< boost::shared_ptr<const Foo> >();
      python::class_<std::vector<boost::shared_ptr<const Foo> > >(
        \"FooList\")
        .def(python::vector_indexing_suite<
             std::vector<boost::shared_ptr<const Foo> >,
             true>() );
  }
};
}
int main ()
{ return 0;
}
" DOES_BOOST_PYTHON_SUPPORT_SHARED_CONST)

  if (DOES_BOOST_PYTHON_SUPPORT_SHARED_CONST)
     add_definitions("-DBOOST_PYTHON_SUPPORT_SHARED_CONST")
  endif()
endif()

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
      message("== Updating Filters.cpp from pains file")
      execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/update_pains.py
        RESULT_VARIABLE rc
      )
      if(NOT rc EQUAL 0)
        message(FATAL_ERROR "Failed updating Filters.cpp from wehi_pains.csv")
      endif(NOT rc EQUAL 0)
endif(PYTHONINTERP_FOUND)

rdkit_library(FilterCatalog
              Filters.cpp
              FilterCatalog.cpp
              FilterCatalogEntry.cpp
              FilterMatchers.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral ${Boost_SERIALIZATION_LIBRARY} )

rdkit_headers(FilterCatalogEntry.h
              FilterCatalog.h
              FilterMatcherBase.h
              FilterMatchers.h
              DEST GraphMol/FilterCatalog)

add_subdirectory(Wrap)

rdkit_test(filterCatalogTest filtercatalogtest.cpp
           LINK_LIBRARIES FilterCatalog 
           Catalogs SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral 
           ${Boost_SERIALIZATION_LIBRARY})
