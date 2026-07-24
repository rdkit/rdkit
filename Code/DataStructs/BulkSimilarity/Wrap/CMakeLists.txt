remove_definitions(-DRDKIT_BULKSIMILARITY_BUILD)
rdkit_python_extension(cBulkSimilarity
                       rdBulkSimilarity.cpp
                       DEST DataStructs
                       LINK_LIBRARIES BulkSimilarity DataStructs)

add_pytest(pyBulkSimilarity
           ${CMAKE_CURRENT_SOURCE_DIR}/testBulkSimilarity.py)
