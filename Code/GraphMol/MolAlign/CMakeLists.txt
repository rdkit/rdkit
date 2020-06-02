
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLALIGN_BUILD)
rdkit_library(MolAlign AlignMolecules.cpp O3AAlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment GraphMol RDGeneral ForceFieldHelpers)

rdkit_headers(AlignMolecules.h O3AAlignMolecules.h DEST GraphMol/MolAlign)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign Descriptors
FileParsers DistGeomHelpers SmilesParse)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()