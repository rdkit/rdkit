rdkit_test(testMMFFForceFieldHelpers testMMFFHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers FileParsers MolTransforms SmilesParse)

rdkit_test(testMMFFForceFieldHelpersMultiThread testMultiThread.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers FileParsers MolTransforms SmilesParse)
