rdkit_library(StructChecker
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp Utilites.cpp Pattern.cpp Stereo.cpp AtomSymbolMatch.cpp
              LINK_LIBRARIES Depictor FileParsers ChemTransforms)

rdkit_headers(StructChecker.h Utilites.h Pattern.h Stereo.h
              DEST GraphMol/StructChecker)

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker ChemTransforms Depictor FileParsers SmilesParse 
                           GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
