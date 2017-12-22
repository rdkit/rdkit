rdkit_test(testUFFForceFieldHelpers testUFFHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           ${RDK_COORDGEN_LIBS} Depictor 
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse
	   SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment
	   GraphMol RDGeometryLib RDGeneral  ${RDKit_THREAD_LIBS})
