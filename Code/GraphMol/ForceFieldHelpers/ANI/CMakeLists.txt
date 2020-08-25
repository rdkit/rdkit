rdkit_catch_test(ANIForceFieldBuilderTest ANIForceFieldHelperTest.cpp catch_main.cpp
    LINK_LIBRARIES 
    DistGeomHelpers ForceFieldHelpers
    FileParsers MolTransforms SmilesParse 
    SubstructMatch MolAlign ANIDescriptors ANIForceField ANIForceFieldHelpers)
