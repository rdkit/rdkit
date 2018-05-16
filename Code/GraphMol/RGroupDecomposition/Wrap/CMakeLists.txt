remove_definitions(-DRDKIT_RGROUPDECOMPOSITION_BUILD)
rdkit_python_extension(rdRGroupDecomposition
                       rdRGroupComposition.cpp
                       DEST Chem 
                       LINK_LIBRARIES
     RGroupDecomposition FMCS ChemTransforms SubstructMatch SmilesParse GraphMol
     RDGeneral RDGeometryLib RDBoost
    ${RDKit_THREAD_LIBS})



#add_pytest(pyRGroupDecomposition
#         ${CMAKE_CURRENT_SOURCE_DIR}/testGroupDecomposition.py)


add_pytest(pyRGroupDecomposition
         ${CMAKE_CURRENT_SOURCE_DIR}/test_rgroups.py)
