include_directories( ${RDKit_ExternalDir}/AvalonTools )
rdkit_python_extension(pyAvalonTools
                       pyAvalonTools.cpp
                       DEST Avalon
                       LINK_LIBRARIES
                       AvalonLib avalon_clib
		       SubstructMatch
		       SmilesParse FileParsers GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
add_pytest(pyAvalonTools
         ${CMAKE_CURRENT_SOURCE_DIR}/testAvalonTools.py)
