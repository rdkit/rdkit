remove_definitions(-DRDKIT_RGROUPDECOMPOSITION_BUILD)
rdkit_python_extension(rdRGroupDecomposition
                       rdRGroupComposition.cpp
                       DEST Chem 
                       LINK_LIBRARIES
     RGroupDecomposition )



#add_pytest(pyRGroupDecomposition
#         ${CMAKE_CURRENT_SOURCE_DIR}/testGroupDecomposition.py)


add_pytest(pyRGroupDecomposition
         ${CMAKE_CURRENT_SOURCE_DIR}/test_rgroups.py)

add_jupytertest(jupyterRGroupDecomposition
         ${CMAKE_CURRENT_SOURCE_DIR}/../test_data rgd_testing.ipynb)