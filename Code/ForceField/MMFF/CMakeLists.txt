if(RDK_TEST_MMFF_COMPLIANCE)
  if((NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_data/MMFF94_reference.log)
    OR (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_data/MMFF94s_reference.log))
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/test_data/MMFF_reference_logs.tgz
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_data
    )
  endif()
rdkit_test(testMMFFForceField testMMFFForceField.cpp 
    LINK_LIBRARIES 
    DistGeomHelpers DistGeometry ForceFieldHelpers ForceField 
    FileParsers MolTransforms SmilesParse SubstructMatch MolAlign
    Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral)
endif(RDK_TEST_MMFF_COMPLIANCE)

