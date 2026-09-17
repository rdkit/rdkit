rdkit_nanobind_extension(rdMetricMatrixCalc rdMetricMatrixCalc.cpp
                       DEST DataManip/Metric
                       LINK_LIBRARIES
                       DataStructs)

add_pytest(pyMatCalc ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMatricCalc.py)
