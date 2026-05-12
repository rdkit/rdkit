remove_definitions(-DRDKIT_MOLCATALOG_BUILD)
rdkit_python_extension(rdMolCatalog rdMolCatalog.cpp
                       DEST Chem
                       LINK_LIBRARIES MolCatalog )

add_pytest(pyMolCatalog ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)



