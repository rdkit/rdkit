if (MSVC)
  add_definitions( "/D_CRT_SECURE_NO_WARNINGS")
endif(MSVC)

rdkit_library(MMPA
              MMPA.cpp
              LINK_LIBRARIES SubstructMatch Depictor FileParsers ChemTransforms GraphMol ${RDKit_THREAD_LIBS} )

rdkit_headers(MMPA.h
              DEST GraphMol/MMPA)

rdkit_test(testMMPA MMPA_UnitTest.cpp LINK_LIBRARIES
  MMPA ChemTransforms Depictor FileParsers SmilesParse 
  GraphMol RDGeneral RDGeometryLib SubstructMatch  ${RDKit_THREAD_LIBS})


add_subdirectory(Wrap)
