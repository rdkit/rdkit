if(RDK_HAS_EIGEN)
rdkit_catch_test(ANIForceFieldCatchTest ANIForceField_catch.cpp  catch_main.cpp
    LINK_LIBRARIES 
    DistGeomHelpers ForceFieldHelpers
    FileParsers MolTransforms SmilesParse 
    SubstructMatch MolAlign ANIDescriptors)
endif(RDK_HAS_EIGEN)