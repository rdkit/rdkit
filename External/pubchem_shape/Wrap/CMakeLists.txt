find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(Python COMPONENTS Interpreter REQUIRED)
set(py3lib "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
find_package(Boost COMPONENTS ${py3lib} REQUIRED CONFIG)


rdkit_python_extension(rdShapeAlign
                       cshapealign.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       PubChemShape )
add_pytest(pyShapeAlign
         ${CMAKE_CURRENT_SOURCE_DIR}/test_rdshapealign.py)
