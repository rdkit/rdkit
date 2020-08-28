
rdkit_library(SimDivPickers
              DistPicker.cpp MaxMinPicker.cpp HierarchicalClusterPicker.cpp
              LINK_LIBRARIES hc RDGeneral)
target_compile_definitions(SimDivPickers PRIVATE RDKIT_SIMDIVPICKERS_BUILD)

rdkit_headers(DistPicker.h LeaderPicker.h 
              HierarchicalClusterPicker.h
              MaxMinPicker.h DEST SimDivPickers)

rdkit_test(testSimDivPickers testPickers.cpp LINK_LIBRARIES SimDivPickers)

rdkit_catch_test(pickersTestsCatch catch_tests.cpp
           LINK_LIBRARIES SimDivPickers DataStructs)


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
