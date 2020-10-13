remove_definitions(-DRDKIT_FORCEFIELDHELPERS_BUILD)
rdkit_test(testUFFForceFieldHelpers testUFFHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers FileParsers MolTransforms SmilesParse)
