rdkit_python_extension(rdShape
                       rdShape.cpp
                       DEST Chem 
                       LINK_LIBRARIES  GraphMol RDGeneral RDGeometryLib DataStructs Shape FileParsers RDGeneral RDGeometryLib RDBoost)

