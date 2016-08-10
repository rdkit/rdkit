rdkit_library(StructChecker
              AugmentedAtomData.cpp 
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp Utilites.cpp
              Pattern.cpp Stereo.cpp AtomSymbolMatch.cpp StripSmallFragments.cpp ReCharge.cpp
              LINK_LIBRARIES Depictor SubstructMatch FileParsers ChemTransforms Descriptors)

rdkit_headers(StructChecker.h Utilites.h Pattern.h Stereo.h StripSmallFragments.h ReCharge.h
              DEST GraphMol/StructChecker)

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker ChemTransforms Depictor Descriptors
                          SubstructMatch FileParsers SmilesParse 
                          GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
