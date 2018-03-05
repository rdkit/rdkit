rdkit_library(SubstructMatch 
              SubstructMatch.cpp SubstructUtils.cpp
              LINK_LIBRARIES GraphMol
                ${RDKit_THREAD_LIBS} )

rdkit_headers(SubstructMatch.h
              SubstructUtils.h DEST GraphMol/Substruct)

rdkit_test(testSubstructMatch test1.cpp LINK_LIBRARIES  FileParsers SmilesParse SubstructMatch
GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS} )

