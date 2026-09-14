if(RDK_USE_BOOST_IOSTREAMS)

rdkit_library(DASHTree
              AtomFeatures.cpp DASHTree.cpp Charges.cpp
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(DASHTree PRIVATE RDKIT_DASHTREE_BUILD)

rdkit_headers(AtomFeatures.h
              DASHTree.h
              DEST GraphMol/DASHTree)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

else(RDK_USE_BOOST_IOSTREAMS)
message("The DASHTree module memory-maps its data files via boost::iostreams, so not building.")
endif(RDK_USE_BOOST_IOSTREAMS)
