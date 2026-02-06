find_package(Python3 COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(py3lib "python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")
find_package(Boost COMPONENTS ${py3lib} REQUIRED CONFIG)


rdkit_python_extension(rdGaussianShape
        rdGaussianShape.cpp
        DEST Chem
        LINK_LIBRARIES GaussianShape)
add_pytest(pyShapeAlign
        ${CMAKE_CURRENT_SOURCE_DIR}/test_gaussian_shape.py)
