remove_definitions(-DRDKIT_DISTGEOMHELPERS_BUILD)
ADD_DEFINITIONS("-DBOOST_PYTHON_MAX_ARITY=32")

rdkit_python_extension(rdDistGeom rdDistGeom.cpp
                       DEST Chem
                       LINK_LIBRARIES
           LINK_LIBRARIES
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse SubstructMatch ForceFieldHelpers ForceField MolAlign Optimizer EigenSolvers Alignment GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)

add_pytest(pyDistGeomHelpers ${CMAKE_CURRENT_SOURCE_DIR}/testDistGeom.py)
