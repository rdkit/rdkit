rdkit_python_extension(rdGaussianShape
        rdGaussianShape.cpp
        DEST Chem
        LINK_LIBRARIES GaussianShape)
add_pytest(pyShapeAlign
        ${CMAKE_CURRENT_SOURCE_DIR}/test_rdgaussian_shape.py)
