rdkit_library(ShapeHelpers SHARED ShapeEncoder.cpp ShapeUtils.cpp
              LINK_LIBRARIES MolTransforms)

rdkit_test(testShapeHelpers testShapeHelpers.cpp
           LINK_LIBRARIES ShapeHelpers FileParsers MolAlign)

add_subdirectory(Wrap)

