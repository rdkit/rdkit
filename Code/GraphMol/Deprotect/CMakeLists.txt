remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_DEPROTECT_BUILD)
rdkit_library(Deprotect
              Deprotect.cpp
              LINK_LIBRARIES ChemReactions 
	      FilterCatalog Descriptors Fingerprints DataStructs Depictor
              FileParsers SubstructMatch ChemTransforms GraphMol ${RDKit_SERIALIZATION_LIBS})


rdkit_headers(Deprotect.h
              DEST GraphMol/Deprotect)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

rdkit_catch_test(deprotectTest deprotectTest.cpp LINK_LIBRARIES Deprotect )

