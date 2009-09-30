rdkit_python_extension(rdDepictor
  rdDepictor.cpp
  DEST Chem 
  LINK_LIBRARIES Depictor GraphMol 
     RDGeometryLib RDGeneral RDBoost )
add_test( pyDepictor python ${CMAKE_CURRENT_SOURCE_DIR}/testDepictor.py )
