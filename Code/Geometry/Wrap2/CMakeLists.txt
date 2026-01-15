remove_definitions(-DRDKIT_RDGEOMETRYLIB_BUILD)


rdkit_nanobind_extension(rdGeometry
    rdGeometry.cpp Point.cpp DEST Geometry
    LINK_LIBRARIES RDGeometryLib rdkit_base)
rdkit_nanobind_extension(rdGeometryTest
    rdGeometryTest.cpp DEST Geometry
    LINK_LIBRARIES RDGeometryLib rdkit_base)


#-------

add_pytest(testPyGeometry ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py)

