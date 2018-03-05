rdkit_python_extension(rdConformerParser rdConformerParser.cpp
                       DEST Chem
			           LINK_LIBRARIES
			           ConformerParser
			           FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral RDBoost
			           ) 

add_pytest(pyConformerParser ${CMAKE_CURRENT_SOURCE_DIR}/testConformerParser.py)
