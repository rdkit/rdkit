include_directories( ${RDKit_ExternalDir}/FreeSASA )
rdkit_python_extension(rdFreeSASA 
                       rdFreeSASA.cpp 
                       DEST Chem
                       LINK_LIBRARIES 
                       FreeSASALib freesasa_clib GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
add_pytest(pyFreeSASA
         ${CMAKE_CURRENT_SOURCE_DIR}/testFreeSASA.py)
