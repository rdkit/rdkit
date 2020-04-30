if(NOT LIB_FUZZING_ENGINE)
  message("LIB_FUZZING_ENGINE is not set, using standalone fuzz runner, no actual fuzzing will be done!")
  add_library(standalone_fuzz_runner STATIC standalone_fuzz_target_runner.cpp)
  set(LIB_FUZZING_ENGINE standalone_fuzz_runner)
endif()

include_directories ( "include" )

set(RDKit_LIBS
           ChemReactions_static FileParsers_static SmilesParse_static Depictor_static
           RDGeometryLib_static RDGeneral_static SubstructMatch_static Subgraphs_static
           MolDraw2D_static GraphMol_static RingDecomposerLib_static DistGeometry_static
           DistGeomHelpers_static MolAlign_static Optimizer_static ForceField_static
           ForceFieldHelpers_static Alignment_static ForceField_static MolTransforms_static
           EigenSolvers_static ForceFieldHelpers_static Alignment_static MolTransforms_static
           EigenSolvers_static
)

macro(fuzzer name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} rdkit_base ${RDKit_LIBS} ${LIB_FUZZING_ENGINE} "/usr/local/lib/libboost_serialization.a")
endmacro()

fuzzer(smiles_string_to_mol_fuzzer)
fuzzer(mol_block_to_mol_fuzzer)
fuzzer(mol_data_stream_to_mol_fuzzer)