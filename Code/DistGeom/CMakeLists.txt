### ISSUE: This library was originally named DistGeom. I've renamed it here
### to avoid a target name clash in ./Wrap

rdkit_library(DistGeometry
              DistGeomUtils.cpp TriangleSmooth.cpp DistViolationContrib.cpp DistViolationContribs.cpp 
              ChiralViolationContrib.cpp ChiralViolationContribs.cpp
              LINK_LIBRARIES EigenSolvers ForceField ForceFieldHelpers)
target_compile_definitions(DistGeometry PRIVATE RDKIT_DISTGEOMETRY_BUILD)

rdkit_headers(BoundsMatrix.h
              ChiralSet.h
              ChiralViolationContrib.h
              ChiralViolationContribs.h
              DistGeomUtils.h
              DistViolationContrib.h
              DistViolationContribs.h
              FourthDimContrib.h
              FourthDimContribs.h
              TriangleSmooth.h DEST DistGeom)

rdkit_test(testDistGeom testDistGeom.cpp 
LINK_LIBRARIES DistGeometry )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()