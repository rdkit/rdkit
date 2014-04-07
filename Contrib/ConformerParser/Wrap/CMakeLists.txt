
set(RDK_BUILD_PYTHON_WRAPPERS ON)
set(RDK_PYTHON_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/rdkit")
set(RDKit_PythonDir ${PROJECT_SOURCE_DIR}/../../rdkit)

rdkit_python_extension(rdConformerParser rdConformerParser.cpp
                       DEST Chem
			           LINK_LIBRARIES
			           ConformerParser
			           ${RDKIT_LIBRARIES}
			           ) 

add_pytest(pyConformerParser ${CMAKE_CURRENT_SOURCE_DIR}/testConformerParser.py)

