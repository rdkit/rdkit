
rdkit_library(EnumerateStereoisomers
      EnumerateStereoisomers.cpp Flippers.cpp
      LINK_LIBRARIES SmilesParse FileParsers GraphMol CIPLabeler DistGeomHelpers DistGeometry)

rdkit_headers(EnumerateStereoisomers.h)

rdkit_catch_test(testEnumerateStereoisomers catch_tests.cpp
      LINK_LIBRARIES EnumerateStereoisomers)
