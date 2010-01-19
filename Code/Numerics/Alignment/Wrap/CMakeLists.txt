rdkit_python_extension(rdAlignment
                       rdAlignment.cpp
                       DEST Numerics 
                       LINK_LIBRARIES Alignment RDGeneral RDGeometryLib RDBoost)

add_test(pyAlignment ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testAlignment.py)

