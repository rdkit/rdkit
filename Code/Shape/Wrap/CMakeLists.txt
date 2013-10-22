rdkit_python_extension(rdShape
                       rdShape.cpp
                       DEST Chem 
                       LINK_LIBRARIES  GraphMol RDGeneral RDGeometryLib DataStructs Shape FileParsers RDGeneral RDGeometryLib RDBoost)

add_pytest(pyShapeReference ${CMAKE_CURRENT_SOURCE_DIR}/test_shape_reference.py)
add_pytest(pyShapeSanity ${CMAKE_CURRENT_SOURCE_DIR}/test_shape_sanity.py)

