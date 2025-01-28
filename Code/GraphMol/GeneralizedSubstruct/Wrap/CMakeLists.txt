rdkit_python_extension(rdGeneralizedSubstruct
                       rdGeneralizedSubstruct.cpp
                       DEST Chem
                       LINK_LIBRARIES GeneralizedSubstruct 
GraphMol ) 

add_pytest(pyGeneralizedSubstruct ${CMAKE_CURRENT_SOURCE_DIR}/testGeneralizedSubstruct.py)

