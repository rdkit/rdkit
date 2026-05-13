rdkit_python_extension(rdMetricMatrixCalc rdMetricMatrixCalc.cpp
                       DEST DataManip/Metric
                       LINK_LIBRARIES
                       DataStructs)

add_pytest(pyMatCalc ${CMAKE_CURRENT_SOURCE_DIR}/testMatricCalc.py)




