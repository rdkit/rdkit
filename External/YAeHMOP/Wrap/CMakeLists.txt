
rdkit_python_extension(rdEHTTools
                       rdEHTTools.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       EHTLib yaehmop_eht 
	                   GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
                       
add_pytest(pyEHTTools
         ${CMAKE_CURRENT_SOURCE_DIR}/testEHTTools.py)
