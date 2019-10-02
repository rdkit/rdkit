file(WRITE ${RDK_PYTHON_BINARY_DIR}/rdkit/RDPaths.py
"import os
# unset so to trigger exceptions and track use: RDBaseDir=os.environ['RDBASE']
RDCodeDir=os.path.dirname(__file__)
# not really hard-coded alternative RDCodeDir=os.path.dirname(__file__)
_share = os.path.join(r'${CMAKE_INSTALL_PREFIX}', r'${RDKit_ShareDir}')
RDDataDir=os.path.join(_share,'Data')
RDDocsDir=os.path.join(_share,'Docs')
RDProjDir=os.path.join(_share,'Projects')
RDContribDir=os.path.join(_share,'Contrib')
")

add_pytest(pythonTestDirML
         ${RDK_PYTHON_BINARY_DIR}/rdkit/ML/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/ML )
add_pytest(pythonTestDirDataStructs
         ${RDK_PYTHON_BINARY_DIR}/rdkit/DataStructs/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/DataStructs )
add_pytest(pythonTestDirDbase
         ${RDK_PYTHON_BINARY_DIR}/rdkit/Dbase/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/Dbase )
add_pytest(pythonTestDirSimDivFilters
         ${RDK_PYTHON_BINARY_DIR}/rdkit/SimDivFilters/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/SimDivFilters )
add_pytest(pythonTestDirVLib
         ${RDK_PYTHON_BINARY_DIR}/rdkit/VLib/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/VLib )
add_pytest(pythonTestDirChem
         ${RDK_PYTHON_BINARY_DIR}/rdkit/Chem/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/Chem --buildType ${CMAKE_BUILD_TYPE})
add_pytest(pythonTestSping
         ${RDK_PYTHON_BINARY_DIR}/rdkit/Chem/test_list.py --testDir ${RDK_PYTHON_BINARY_DIR}/rdkit/sping )

