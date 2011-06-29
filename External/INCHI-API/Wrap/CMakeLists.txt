rdkit_python_extension(rdinchi
                       pyInchi.cpp
		       DEST Chem
                       LINK_LIBRARIES 
                       RDInchiLib Inchi SubstructMatch Depictor GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
