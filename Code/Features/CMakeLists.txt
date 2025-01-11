rdkit_headers(Feature.h
              DEST Features)

rdkit_test(testBaseFeatures testFeatures.cpp 
LINK_LIBRARIES RDGeneral RDGeometryLib )
