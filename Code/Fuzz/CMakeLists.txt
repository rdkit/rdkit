#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(RDBASE $ENV{RDBASE})
set(LIB_FUZZING_ENGINE $ENV{LIB_FUZZING_ENGINE})
if(NOT LIB_FUZZING_ENGINE)
  add_library(standalone_fuzz_runner STATIC standalone_fuzz_target_runner.cpp)
  set(LIB_FUZZING_ENGINE standalone_fuzz_runner)
endif()

# specify where CMake can find the RDKit libraries
#include_directories ( ${RDBASE}/Code )
#link_directories ( ${RDBASE}/lib )
#include_directories(${RDKIT_INCLUDE_DIR})


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
  target_link_libraries(${name} rdkit_base ${RDKit_LIBS} ${LIB_FUZZING_ENGINE})
endmacro()

fuzzer(smile_fuzz)