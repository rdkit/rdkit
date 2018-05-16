remove_definitions(-DRDKIT_SIMDIVPICKERS_BUILD)
rdkit_python_extension(rdSimDivPickers 
                       MaxMinPicker.cpp HierarchicalClusterPicker.cpp 
                       rdSimDivPickers.cpp 
                       DEST SimDivFilters
                       LINK_LIBRARIES SimDivPickers 
hc DataStructs
RDGeneral RDBoost)

add_pytest(pySimDivPickers ${CMAKE_CURRENT_SOURCE_DIR}/testPickers.py)




