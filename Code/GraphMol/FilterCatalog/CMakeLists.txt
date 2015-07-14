set(T_LIBS ${Boost_LIBRARIES})
find_package(Boost 1.39.0 COMPONENTS serialization)
if (Boost_SERIALIZATION_LIBRARY)
  set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
  message("== Using boost serialization library ${Boost_SERIALIZATION_LIBRARY}")
  ADD_DEFINITIONS("-DRDK_USE_BOOST_SERIALIZATION")
else()
  message("== Making FilterCatalog without boost Serialization support")
endif()

rdkit_library(FilterCatalog
              Filters.cpp
              FilterCatalog.cpp
              FilterCatalogEntry.cpp
              FilterMatchers.cpp
              LINK_LIBRARIES ${Boost_SERIALIZATION_LIBRARY} Subgraphs SubstructMatch SmilesParse GraphMol RDGeometryLib Catalogs RDGeneral )

rdkit_headers(FilterCatalogEntry.h
              FilterCatalog.h
              FilterMatcherBase.h
              FilterMatchers.h
              DEST GraphMol/FilterCatalog)

add_subdirectory(Wrap)

rdkit_test(filterCatalogTest filtercatalogtest.cpp
           LINK_LIBRARIES ${Boost_SERIALIZATION_LIBRARY} FilterCatalog SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral)
