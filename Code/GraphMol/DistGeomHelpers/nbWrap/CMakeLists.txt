remove_definitions(-DRDKIT_DISTGEOMHELPERS_BUILD)

rdkit_nanobind_extension(rdDistGeom rdDistGeom.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         DistGeomHelpers FileParsers MolTransforms SmilesParse SubstructMatch MolAlign EigenSolvers)

add_pytest(pyDistGeomHelpers ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testDistGeom.py)
