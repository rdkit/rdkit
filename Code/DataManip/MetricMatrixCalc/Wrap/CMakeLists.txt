rdkit_python_extension(rdMetricMatrixCalc rdMetricMatrixCalc.cpp
                       DEST DataManip/Metric
                       LINK_LIBRARIES
                       RDGeneral DataStructs RDBoost)

add_test(pyMatCalc ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testMatricCalc.py)




