rdkit_python_extension(rdAlignment
                       rdAlignment.cpp
                       DEST Numerics 
                       LINK_LIBRARIES Alignment RDBoost)

add_test(pyAlignment python ${CMAKE_CURRENT_SOURCE_DIR}/testAlignment.py)

