rdkit_python_extension(DistGeom DistGeom.cpp
                       DEST DistanceGeometry
LINK_LIBRARIES DistGeometry EigenSolvers ForceField Optimizer RDGeneral RDGeometryLib RDBoost )


add_test(pyDistGeom python ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)




