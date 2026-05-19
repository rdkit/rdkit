remove_definitions(-DRDKIT_RGROUPDECOMPOSITION_BUILD)
rdkit_nanobind_extension(rdRGroupDecomposition
                         rdRGroupComposition.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         RGroupDecomposition)

add_pytest(pyRGroupDecomposition
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/test_rgroups.py)
