
add_definitions(-DRDKIT_SIMDIVPICKERS_BUILD)
rdkit_library(SimDivPickers
              DistPicker.cpp MaxMinPicker.cpp HierarchicalClusterPicker.cpp
              LINK_LIBRARIES hc RDGeneral)

rdkit_headers(DistPicker.h
              HierarchicalClusterPicker.h
              MaxMinPicker.h DEST SimDivPickers)

rdkit_test(testSimDivPickers testPickers.cpp LINK_LIBRARIES SimDivPickers RDGeneral)

add_subdirectory(Wrap)
