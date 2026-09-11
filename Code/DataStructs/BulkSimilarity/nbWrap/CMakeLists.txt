remove_definitions(-DRDKIT_BULKSIMILARITY_BUILD)
rdkit_nanobind_extension(cBulkSimilarity
                         rdBulkSimilarity.cpp
                         DEST DataStructs
                         LINK_LIBRARIES
                         BulkSimilarity DataStructs)

# The Boost.Python and nanobind extensions expose the same API, so they
# share a single test file.
add_pytest(pyBulkSimilarity
           ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testBulkSimilarity.py)
