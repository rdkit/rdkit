rdkit_python_extension(rdSimDivPickers 
                       MaxMinPicker.cpp HierarchicalClusterPicker.cpp 
                       rdSimDivPickers.cpp 
                       DEST SimDivFilters
                       LINK_LIBRARIES SimDivPickers 
hc DataStructs
RDGeneral RDBoost)





