if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making SubstructLibrary without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

rdkit_library(TautomerQuery 
              TautomerQuery.cpp
              LINK_LIBRARIES GraphMol RDGeneral MolStandardize ${RDKit_SERIALIZATION_LIBS})
target_compile_definitions(TautomerQuery PRIVATE RDKIT_TAUTOMERQUERY_BUILD)

rdkit_headers(TautomerQuery.h
              DEST GraphMol/TautomerQuery)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()

rdkit_catch_test(tautomerQueryTestCatch catch_tests.cpp LINK_LIBRARIES TautomerQuery )
