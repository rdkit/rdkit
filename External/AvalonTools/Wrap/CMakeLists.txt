include_directories( ${RDKit_ExternalDir}/AvalonTools )
rdkit_python_extension(pyAvalonTools 
                       pyAvalonTools.cpp 
                       DEST Avalon
                       LINK_LIBRARIES 
                       AvalonLib avalon_clib SmilesParse FileParsers GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
add_test(pyAvalonTools ${PYTHON_EXECUTABLE} 
         ${CMAKE_CURRENT_SOURCE_DIR}/testAvalonTools.py)
