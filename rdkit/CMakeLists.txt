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

add_pythonpytest(pythonSourceTests ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory(Chem)
