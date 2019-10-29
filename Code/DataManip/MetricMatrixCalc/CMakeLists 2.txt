### Template library
rdkit_headers(MetricFuncs.h
              MetricMatrixCalc.h DEST DataManip/MetricMatrixCalc)

rdkit_test(testMatCalc testMatCalc.cpp LINK_LIBRARIES RDGeneral)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()


