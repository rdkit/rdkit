### ISSUE: This library was originally named DistGeom. I've renamed it here
### to avoid a target name clash in ./Wrap
rdkit_library(DistGeometry
              DistGeomUtils.cpp TriangleSmooth.cpp DistViolationContrib.cpp 
              ChiralViolationContrib.cpp
              LINK_LIBRARIES EigenSolvers ForceField)

rdkit_headers(BoundsMatrix.h
              ChiralSet.h
              ChiralViolationContrib.h
              DistGeomUtils.h
              DistViolationContrib.h
              FourthDimContrib.h
              TriangleSmooth.h DEST DistGeom)

rdkit_test(testDistGeom testDistGeom.cpp 
LINK_LIBRARIES DistGeometry EigenSolvers ForceField Optimizer RDGeneral RDGeometryLib )

add_subdirectory(Wrap)
