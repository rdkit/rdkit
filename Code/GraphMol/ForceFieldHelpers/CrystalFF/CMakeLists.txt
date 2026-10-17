rdkit_test(testCrystalForceField testCrystalFF.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers FileParsers MolTransforms SmilesParse )

rdkit_catch_test(
    crystalForceFieldCatch
        catch_tests.cpp
    LINK_LIBRARIES
        DistGeomHelpers
        ForceFieldHelpers
        MolTransforms
        SmilesParse
) 