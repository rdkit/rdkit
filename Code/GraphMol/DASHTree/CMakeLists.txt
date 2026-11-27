if(RDK_USE_BOOST_IOSTREAMS AND NOT RDK_BIG_ENDIAN)

rdkit_library(DASHTree
              AtomFeatures.cpp DASHTree.cpp Charges.cpp
              LINK_LIBRARIES GraphMol RDGeneral)
target_compile_definitions(DASHTree PRIVATE RDKIT_DASHTREE_BUILD)

rdkit_headers(DASHTree.h DEST GraphMol/DASHTree)

# the tree the tests and GetDASHTree() use by default, too large for the repo
set(dashTree "${RDKit_DataDir}/DASHTree/initial_test.dash")
set(dashTreeMD5 "3dc6a5c2f221e7a3be2edea1b38518b3")
if(EXISTS "${dashTree}")
  file(MD5 "${dashTree}" foundMD5)
endif()
if(NOT "${foundMD5}" STREQUAL "${dashTreeMD5}")
  downloadAndCheckMD5("https://github.com/MTLehner/dash-trees/raw/v1/initial_test.dash"
                      "${dashTree}" "${dashTreeMD5}")
endif()

rdkit_catch_test(dashTreeCatchTest catch_tests.cpp
                 LINK_LIBRARIES DASHTree SmilesParse)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  add_subdirectory(nbWrap)
endif()

else()
message("The DASHTree module needs boost::iostreams and a little-endian host, so it is not built.")
endif()
