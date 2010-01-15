add_test(pythonTestDirML
         python ${CMAKE_CURRENT_SOURCE_DIR}/ML/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/ML )
add_test(pythonTestDirDataStructs
         python ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs )
add_test(pythonTestDirDbase
         python ${CMAKE_CURRENT_SOURCE_DIR}/Dbase/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Dbase )
add_test(pythonTestDirSimDivFilters
         python ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters )
add_test(pythonTestDirVLib
         python ${CMAKE_CURRENT_SOURCE_DIR}/VLib/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/VLib )
add_test(pythonTestDirChem
         python ${CMAKE_CURRENT_SOURCE_DIR}/Chem/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Chem )

