rdkit_python_extension(rdNormalizedDescriptors
                       rdNormalizedDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       NormalizedDescriptors)

add_pytest(pyNormalizedDescriptors
         ${CMAKE_CURRENT_SOURCE_DIR}/testNormalizedDescriptors.py)
