remove_definitions(-DRDKIT_SYNTHONSPACESEARCH_BUILD)

rdkit_python_extension(rdSynthonSpaceSearch
        rdSynthonSpaceSearch.cpp
        DEST Chem
        LINK_LIBRARIES SynthonSpaceSearch)

#add_pytest(pyHyperspaceSearch ${CMAKE_CURRENT_SOURCE_DIR}/testHyperspaceSearch.py)
