
rdkit_library(MonomerMol
    MonomerMol.cpp toAtomistic.cpp toMonomeric.cpp MonomerLibrary.cpp
    LINK_LIBRARIES GenericGroups Depictor SmilesParse ChemTransforms GraphMol)
target_compile_definitions(FileParsers PRIVATE RDKIT_FILEPARSERS_BUILD)

rdkit_headers(MonomerMol.h Conversions.h
    DEST GraphMol/MonomerMol)

rdkit_catch_test(testMonomerMol testMonomerMol.cpp
    LINK_LIBRARIES MonomerMol FileParsers SubstructMatch GraphMol)

