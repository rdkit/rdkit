remove_definitions(-DRDKIT_TAUTOMERQUERY_BUILD)
rdkit_python_extension(rdTautomerQuery
                       rdTautomerQuery.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       TautomerQuery )

add_pytest(pyTautomerQuery ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
