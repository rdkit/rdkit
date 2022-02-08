file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/RDPaths.py 
"import os
# unset so to trigger exceptions and track use: RDBaseDir=os.environ['RDBASE']
RDCodeDir=os.path.join(r'${PYTHON_INSTDIR}','rdkit')
# not really hard-coded alternative RDCodeDir=os.path.dirname(__file__)
_share = os.path.join(r'${CMAKE_INSTALL_PREFIX}', r'${RDKit_ShareDir}')
RDDataDir=os.path.join(_share,'Data')
RDDocsDir=os.path.join(_share,'Docs')
RDProjDir=os.path.join(_share,'Projects')
RDContribDir=os.path.join(_share,'Contrib')
")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RDPaths.py 
  DESTINATION ${RDKit_PythonDir}
  COMPONENT python)

add_pytest(pythonTestDirRoot
         ${CMAKE_CURRENT_SOURCE_DIR}/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR} )
add_pytest(pythonTestDirML
         ${CMAKE_CURRENT_SOURCE_DIR}/ML/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/ML )
add_pytest(pythonTestDirDataStructs
         ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/DataStructs )
add_pytest(pythonTestDirDbase
         ${CMAKE_CURRENT_SOURCE_DIR}/Dbase/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/Dbase )
add_pytest(pythonTestDirSimDivFilters
         ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/SimDivFilters )
add_pytest(pythonTestDirVLib
         ${CMAKE_CURRENT_SOURCE_DIR}/VLib/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/VLib )
add_pytest(pythonTestSping
         ${CMAKE_CURRENT_SOURCE_DIR}/Chem/test_list.py --testDir ${CMAKE_CURRENT_SOURCE_DIR}/sping )

add_subdirectory(Chem)
