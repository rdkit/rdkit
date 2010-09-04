rdkit_library(SimDivPickers
              DistPicker.cpp MaxMinPicker.cpp HierarchicalClusterPicker.cpp
              LINK_LIBRARIES hc RDGeneral)

rdkit_headers(DistPicker.h
              HierarchicalClusterPicker.h
              MaxMinPicker.h DEST SimDivPickers)

add_subdirectory(Wrap)


