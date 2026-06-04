rdkit_nanobind_extension(rdConformerParser rdConformerParser.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         ConformerParser
                         FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral)

add_pytest(pyConformerParser ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testConformerParser.py)
