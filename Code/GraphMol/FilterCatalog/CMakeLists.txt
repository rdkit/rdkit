if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    ADD_DEFINITIONS("-DRDK_USE_BOOST_SERIALIZATION")
else()
      message("== Making FilterCatalog without boost Serialization support")
endif()

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
  if(${PYTHON_VERSION_STRING} VERSION_GREATER "2.6" )
        message("== Updating Filters.cpp from pains file")
        execute_process(
          COMMAND
          ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/update_pains.py
          RESULT_VARIABLE rc
        )
        if(NOT rc EQUAL 0)
          message(FATAL_ERROR "Failed updating Filters.cpp from wehi_pains.csv")
        endif(NOT rc EQUAL 0)
  endif()
endif()

rdkit_library(FilterCatalog
              Filters.cpp
              FilterCatalog.cpp
              FilterCatalogEntry.cpp
              FilterMatchers.cpp
              FunctionalGroupHierarchy.cpp
              LINK_LIBRARIES Subgraphs SubstructMatch SmilesParse
              GraphMol RDGeometryLib Catalogs RDGeneral ${Boost_SERIALIZATION_LIBRARY} )

rdkit_headers(FilterCatalogEntry.h
              FilterCatalog.h
              FilterMatcherBase.h
              FilterMatchers.h
              FunctionalGroupHierarchy.h
              DEST GraphMol/FilterCatalog)

add_subdirectory(Wrap)

rdkit_test(filterCatalogTest filtercatalogtest.cpp
           LINK_LIBRARIES FilterCatalog
           Catalogs SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral
           ${Boost_SERIALIZATION_LIBRARY})
