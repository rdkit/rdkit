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

if(NOT RDK_INSTALL_INTREE AND RDK_INSTALL_PYTHON_TESTS)
  set(RDKit_PythonSourceTestsDir ${RDKit_PythonDir})
else()
  set(RDKit_PythonSourceTestsDir ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_pythonpytest(pythonSourceTests ${RDKit_PythonSourceTestsDir} )


add_subdirectory(Chem)
