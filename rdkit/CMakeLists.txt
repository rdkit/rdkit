add_test(pythonTestDirML
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ML/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/ML )
add_test(pythonTestDirDataStructs
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs )
add_test(pythonTestDirDbase
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Dbase/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Dbase )
add_test(pythonTestDirSimDivFilters
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters )
add_test(pythonTestDirVLib
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/VLib/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/VLib )
add_test(pythonTestDirChem
         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Chem/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Chem )

