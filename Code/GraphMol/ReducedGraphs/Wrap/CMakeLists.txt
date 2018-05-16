remove_definitions(-DRDKIT_REDUCEDGRAPHS_BUILD)
rdkit_python_extension(rdReducedGraphs rdReducedGraphs.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES ReducedGraphs SubstructMatch SmilesParse GraphMol DataStructs RDGeneral RDGeometryLib RDBoost )

add_pytest(pyReducedGraphs ${CMAKE_CURRENT_SOURCE_DIR}/testReducedGraphs.py)





