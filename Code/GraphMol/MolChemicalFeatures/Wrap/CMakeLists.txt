rdkit_python_extension(rdMolChemicalFeatures
                       rdMolChemicalFeatures.cpp MolChemicalFeature.cpp 
                       MolChemicalFeatureFactory.cpp ChemicalFeatureUtils.cpp
		       DEST Chem
		       LINK_LIBRARIES MolChemicalFeatures 
Catalogs
SmilesParse SubstructMatch DataStructs
GraphMol RDGeometryLib RDGeneral
RDBoost)

add_pytest(pyChemicalFeatures
         ${CMAKE_CURRENT_SOURCE_DIR}/testChemicalFeatures.py)




