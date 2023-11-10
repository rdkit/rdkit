
if(RDK_BUILD_YAEHMOP_SUPPORT)
set(ehtlib EHTLib)
endif()
rdkit_library(DetermineBonds DetermineBonds.cpp
              LINK_LIBRARIES ${ehtlib} FileParsers GraphMol RDGeneral)

target_compile_definitions(DetermineBonds PRIVATE RDKIT_DETERMINEBONDS_BUILD)

rdkit_headers(DetermineBonds.h DEST GraphMol/DetermineBonds)

rdkit_catch_test(determineBondsCatchTest catch_tests.cpp 
           LINK_LIBRARIES DetermineBonds FileParsers SmilesParse DistGeomHelpers ForceFieldHelpers)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()           