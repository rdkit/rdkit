rdkit_library(GaussianShape
        GaussianShape.cpp ShapeInput.cpp SingleConformerAlignment.cpp
        SHARED LINK_LIBRARIES SmilesParse SubstructMatch MolTransforms SimDivPickers)
target_compile_definitions(GaussianShape PRIVATE RDKIT_GAUSSIANSHAPE_BUILD)

rdkit_headers(GaussianShape.h ShapeInput.h ShapeOverlayOptions.h)

rdkit_catch_test(testGaussianShape catch_tests.cpp LINK_LIBRARIES GaussianShape
        FileParsers MolAlign MolTransforms DistGeomHelpers DistGeometry)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()
