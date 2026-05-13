remove_definitions(-DRDKIT_ALIGNMENT_BUILD)
rdkit_python_extension(rdAlignment
                       rdAlignment.cpp
                       DEST Numerics 
                       LINK_LIBRARIES Alignment)

add_pytest(pyAlignment ${CMAKE_CURRENT_SOURCE_DIR}/testAlignment.py)

