remove_definitions(-DRDKIT_GAUSSIANSHAPE_BUILD)
rdkit_nanobind_extension(rdGaussianShape
        rdGaussianShape.cpp
        DEST Chem
        LINK_LIBRARIES GaussianShape)
add_pytest(pyShapeAlign
        ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/test_rdgaussian_shape.py)
