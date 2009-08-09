rdkit_python_extension(rdchemtransforms
  rdchemtransforms.cpp
  DEST Chem 
  LINK_LIBRARIES ChemTransforms SubstructMatch GraphMol 
     RDGeometry DataStructs RDGeneral RDBoost )
add_test( pyChemTransforms python ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py )

