remove_definitions(-DRDKIT_MOLCHEMICALFEATURES_BUILD)
rdkit_nanobind_extension(rdMolChemicalFeatures
                         rdMolChemicalFeatures.cpp MolChemicalFeature.cpp
                         MolChemicalFeatureFactory.cpp ChemicalFeatureUtils.cpp
                         DEST Chem
                         LINK_LIBRARIES MolChemicalFeatures)

add_pytest(pyChemicalFeatures
           ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testChemicalFeatures.py)
