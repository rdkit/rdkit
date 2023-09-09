
rdkit_library(GeneralizedSubstruct 
              XQMol.cpp TextIO.cpp
              LINK_LIBRARIES MolEnumerator TautomerQuery SubstructMatch SmilesParse GraphMol) 
target_compile_definitions(GeneralizedSubstruct PRIVATE RDKIT_GENERALIZEDSUBSTRUCT_BUILD)

rdkit_headers(XQMol.h DEST GraphMol/GeneralizedSubstruct)

rdkit_catch_test(generalizedSubstructCatch catch_tests.cpp LINK_LIBRARIES GenericGroups FileParsers GeneralizedSubstruct)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
