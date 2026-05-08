remove_definitions(-DRDKIT_PARTIALCHARGES_BUILD)
rdkit_python_extension(rdPartialCharges
                       rdPartialCharges.cpp
                       DEST Chem 
                       LINK_LIBRARIES PartialCharges ) 

add_pytest(pyPartialCharges ${CMAKE_CURRENT_SOURCE_DIR}/testPartialCharges.py)
