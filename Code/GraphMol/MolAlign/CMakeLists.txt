rdkit_library(MolAlign AlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign FileParsers DistGeomHelpers)

add_subdirectory(Wrap)




