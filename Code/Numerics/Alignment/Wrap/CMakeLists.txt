rdkit_python_extension(rdAlignment
                       rdAlignment.cpp
                       DEST Numerics 
                       LINK_LIBRARIES Alignment RDGeneral RDGeometryLib RDBoost)

add_pytest(pyAlignment ${CMAKE_CURRENT_SOURCE_DIR}/testAlignment.py)

