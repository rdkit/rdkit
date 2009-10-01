rdkit_python_extension(rdShapeHelpers rdShapeHelpers.cpp
                       DEST Chem
		       LINK_LIBRARIES
                       ShapeHelpers 
MolTransforms EigenSolvers
DataStructs GraphMol RDGeometryLib RDGeneral
RDBoost)

add_test(pyShapeHelpers python
         ${CMAKE_CURRENT_SOURCE_DIR}/testShapeHelpers.py)


