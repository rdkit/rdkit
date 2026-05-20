remove_definitions(-DRDKIT_DEPICTOR_BUILD)
rdkit_python_extension(rdDepictor
                       rdDepictor.cpp
                       DEST Chem
                       LINK_LIBRARIES
            Depictor )

add_pytest(pyDepictor ${CMAKE_CURRENT_SOURCE_DIR}/testDepictor.py)
