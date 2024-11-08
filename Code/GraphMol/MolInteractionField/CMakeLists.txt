rdkit_library(RDMIFLib MIFDescriptors.cpp LINK_LIBRARIES ForceFieldHelpers ForceField DataStructs RDGeneral GraphMol SmilesParse)

rdkit_headers(MIFDescriptors.h DEST GraphMol/MIF)

rdkit_test(testMIF testMIF.cpp LINK_LIBRARIES RDMIFLib DistGeomHelpers DistGeometry ForceFieldHelpers ForceField
 PartialCharges SubstructMatch FileParsers SmilesParse
MolAlign MolTransforms
GraphMol
Optimizer EigenSolvers
Alignment
DataStructs RDGeometryLib RDGeneral
${RDKit_THREAD_LIBS} )

add_subdirectory(Wrap)

