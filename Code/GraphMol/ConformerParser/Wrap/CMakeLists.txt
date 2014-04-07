rdkit_python_extension(rdConfParser rdConfParser.cpp
                       DEST Chem
                       LINK_LIBRARIES
			           LINK_LIBRARIES
			           ConformerParser FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral RDBoost) 

add_pytest(pyConfParser ${CMAKE_CURRENT_SOURCE_DIR}/testConfParser.py)

