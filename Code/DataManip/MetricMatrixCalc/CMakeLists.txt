### Template library
rdkit_headers(MetricFuncs.h
              MetricMatrixCalc.h DEST DataManip/MetricMatrixCalc)

rdkit_test(testMatCalc testMatCalc.cpp LINK_LIBRARIES RDGeneral)

add_subdirectory(Wrap)



