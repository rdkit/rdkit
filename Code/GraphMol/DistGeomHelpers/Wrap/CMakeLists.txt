rdkit_python_extension(rdDistGeom rdDistGeom.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       DistGeomHelpers RDBoost) 

add_test(pyDistGeom python ${CMAKE_CURRENT_SOURCE_DIR}/testDistGeom.py)

