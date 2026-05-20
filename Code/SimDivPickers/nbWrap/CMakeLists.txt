remove_definitions(-DRDKIT_SIMDIVPICKERS_BUILD)
rdkit_nanobind_extension(rdSimDivPickers
                       MaxMinPicker.cpp LeaderPicker.cpp HierarchicalClusterPicker.cpp
                       rdSimDivPickers.cpp
                       DEST SimDivFilters
                       LINK_LIBRARIES SimDivPickers DataStructs)

add_pytest(pySimDivPickersMaxMin ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMaxMin.py)
add_pytest(pySimDivPickers ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testPickers.py)
