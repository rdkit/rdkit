if (MSVC)
  add_definitions( "/D_CRT_SECURE_NO_WARNINGS")
endif(MSVC)


rdkit_library(MMPA
              MMPA.cpp
              LINK_LIBRARIES SubstructMatch Depictor FileParsers ChemTransforms GraphMol RDGeneral )
target_compile_definitions(MMPA PRIVATE RDKIT_MMPA_BUILD)

rdkit_headers(MMPA.h
              DEST GraphMol/MMPA)

rdkit_test(testMMPA MMPA_UnitTest.cpp LINK_LIBRARIES
  MMPA SmilesParse )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
