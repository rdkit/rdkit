add_test(pyTestDirML
         python ${CMAKE_CURRENT_SOURCE_DIR}/ML/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/ML )
add_test(pyTestDirDataStructs
         python ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs )
add_test(pyTestDirDbase
         python ${CMAKE_CURRENT_SOURCE_DIR}/Dbase/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Dbase )
add_test(pyTestDirSimDivFilters
         python ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters )
add_test(pyTestDirVLib
         python ${CMAKE_CURRENT_SOURCE_DIR}/VLib/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/VLib )
add_test(pyTestDirChem
         python ${CMAKE_CURRENT_SOURCE_DIR}/Chem/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Chem )

