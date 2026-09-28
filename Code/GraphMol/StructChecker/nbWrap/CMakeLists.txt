remove_definitions(-DRDKIT_STRUCTCHECKER_BUILD)
rdkit_nanobind_extension(rdStructChecker
        structchecker.cpp
        DEST Chem
        LINK_LIBRARIES StructChecker )
add_pytest(pyStructChecker ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
