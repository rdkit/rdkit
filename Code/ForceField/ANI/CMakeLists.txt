if(RDK_BUILD_ANI)
rdkit_catch_test(ANIForceFieldCatchTest ANIForceField_catch.cpp  catch_main.cpp
    LINK_LIBRARIES 
    DistGeomHelpers ForceFieldHelpers
    FileParsers MolTransforms SmilesParse 
    SubstructMatch MolAlign ANIDescriptors ANIForceField)
endif(RDK_BUILD_ANI)
