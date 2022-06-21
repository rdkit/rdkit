
rdkit_library(MolAlign AlignMolecules.cpp 
              LINK_LIBRARIES MolTransforms SmilesParse SubstructMatch Alignment GraphMol RDGeneral)
target_compile_definitions(MolAlign PRIVATE RDKIT_MOLALIGN_BUILD)

rdkit_library(O3AAlign O3AAlignMolecules.cpp
              LINK_LIBRARIES MolAlign MolTransforms SubstructMatch Alignment GraphMol RDGeneral ForceFieldHelpers)
target_compile_definitions(O3AAlign PRIVATE RDKIT_MOLALIGN_BUILD)

rdkit_headers(AlignMolecules.h O3AAlignMolecules.h DEST GraphMol/MolAlign)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign Descriptors
FileParsers DistGeomHelpers SmilesParse)

rdkit_test(testO3AAlign testO3AAlign.cpp
           LINK_LIBRARIES O3AAlign Descriptors
FileParsers DistGeomHelpers SmilesParse)

rdkit_catch_test(molAlignCatchTest catch_tests.cpp 
           LINK_LIBRARIES MolAlign SmilesParse FileParsers )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
