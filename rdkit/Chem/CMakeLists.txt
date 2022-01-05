add_pytest(pythonTestDirChem
         ${CMAKE_CURRENT_SOURCE_DIR}/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR} --buildType ${CMAKE_BUILD_TYPE})

add_jupytertest(github4823 ${CMAKE_CURRENT_SOURCE_DIR}/nbtests github4823.ipynb)