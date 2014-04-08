
set(RDK_BUILD_PYTHON_WRAPPERS ON)
set(RDKit_PythonDir ${PROJECT_SOURCE_DIR}/../../rdkit)

rdkit_python_extension(rdConformerParser rdConformerParser.cpp
                       DEST Chem
			           LINK_LIBRARIES
			           ConformerParser
			           ${RDKIT_LIBRARIES}
			           ) 

#add_pytest(pyConformerParser ${CMAKE_CURRENT_SOURCE_DIR}/testConformerParser.py)
add_test(pyConformerParser ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testConformerParser.py)

