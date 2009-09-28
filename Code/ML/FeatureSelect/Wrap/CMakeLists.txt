include_directories(${CMAKE_SOURCE_DIR}/External/cmim-1.0)

rdkit_python_extension(rdFeatSelect rdFeatSelect.cpp
                       DEST ML/FeatureSelect
                       LINK_LIBRARIES
                       fastentropy DataStructs RDBoost)

add_test(pyCMIM python ${CMAKE_CURRENT_SOURCE_DIR}/testCMIM.py)
