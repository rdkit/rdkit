rdkit_library(SimDivPickers
              DistPicker.cpp MaxMinPicker.cpp HierarchicalClusterPicker.cpp
              LINK_LIBRARIES hc RDGeneral)

add_subdirectory(Wrap)


