remove_definitions(-DRDKIT_DISTGEOMETRY_BUILD)
rdkit_nanobind_extension(DistGeom DistGeom.cpp
                         DEST DistanceGeometry
LINK_LIBRARIES DistGeometry )

add_pytest(pyDistGeom ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
