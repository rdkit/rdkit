remove_definitions(-DRDKIT_MMPA_BUILD)
rdkit_python_extension(rdMMPA
                       rdMMPA.cpp
                       DEST Chem
                       LINK_LIBRARIES MMPA )

add_pytest(pyMMPA ${CMAKE_CURRENT_SOURCE_DIR}/testMMPA.py)
