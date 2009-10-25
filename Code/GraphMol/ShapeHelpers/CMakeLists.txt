rdkit_library(ShapeHelpers ShapeEncoder.cpp ShapeUtils.cpp
              LINK_LIBRARIES MolTransforms)

rdkit_test(testShapeHelpers testShapeHelpers.cpp
           LINK_LIBRARIES ShapeHelpers FileParsers MolAlign
SmilesParse FileParsers MolTransforms SubstructMatch EigenSolvers Alignment DataStructs
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)

