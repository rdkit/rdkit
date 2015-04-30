rdkit_library(PartialCharges
              GasteigerCharges.cpp GasteigerParams.cpp
              LINK_LIBRARIES GraphMol)

rdkit_headers(GasteigerCharges.h
              GasteigerParams.h DEST GraphMol/PartialCharges)

rdkit_test(testPartialCharges test1.cpp LINK_LIBRARIES  PartialCharges
SmilesParse GraphMol RDGeometryLib RDGeneral DataStructs 
${RDKit_THREAD_LIBS} )

add_subdirectory(Wrap)
