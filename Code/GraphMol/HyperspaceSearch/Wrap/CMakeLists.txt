remove_definitions(-DRDKIT_HYPERSPACESEARCH_BUILD)

rdkit_python_extension(rdHyperspaceSearch
        rdHyperspaceSearch.cpp
        DEST Chem
        LINK_LIBRARIES HyperspaceSearch)

#add_pytest(pyHyperspaceSearch ${CMAKE_CURRENT_SOURCE_DIR}/testHyperspaceSearch.py)
