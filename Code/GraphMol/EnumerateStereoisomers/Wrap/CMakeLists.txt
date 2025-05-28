remove_definitions(-DRDKIT_ENUMERATESTEREOISOMERS_BUILD)

rdkit_python_extension(rdEnumerateStereoisomers
        rdEnumerateStereoisomers.cpp
        DEST Chem
        LINK_LIBRARIES EnumerateStereoisomers)

add_pytest(pyNewEnumerateStereoisomers ${CMAKE_CURRENT_SOURCE_DIR}/testEnumerateStereoisomers.py)
