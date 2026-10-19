if(RDK_USE_BOOST_IOSTREAMS AND NOT RDK_BIG_ENDIAN)

rdkit_library(DASHTree
              AtomFeatures.cpp DASHTree.cpp Charges.cpp
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(DASHTree PRIVATE RDKIT_DASHTREE_BUILD)

rdkit_headers(DASHTree.h DEST GraphMol/DASHTree)

rdkit_catch_test(dashTreeCatchTest catch_tests.cpp
                 LINK_LIBRARIES DASHTree SmilesParse)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

else()
message("The DASHTree module needs boost::iostreams and a little-endian host, so not building.")
endif()
