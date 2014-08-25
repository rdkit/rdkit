rdkit_test(FMCStest testFMCS.cpp LINK_LIBRARIES
FMCS FileParsers SmilesParse GraphMol SubstructMatch DataStructs RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
