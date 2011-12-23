if(RDK_BUILD_AVALON_SUPPORT)
set(AVALON_SRC_PATH ${AVALONTOOLS_DIR}/common)

add_definitions(-DBUILD_AVALON_SUPPORT)

set(avalon_clib_srcs ${AVALON_SRC_PATH}/layout.c 
	${AVALON_SRC_PATH}/symboltable.c 
	${AVALON_SRC_PATH}/patclean.c 
	${AVALON_SRC_PATH}/utilities.c
	${AVALON_SRC_PATH}/symbol_lists.c
	${AVALON_SRC_PATH}/stereo.c
	${AVALON_SRC_PATH}/set.c
	${AVALON_SRC_PATH}/perceive.c 
	${AVALON_SRC_PATH}/local.c
	${AVALON_SRC_PATH}/graph.c
	${AVALON_SRC_PATH}/geometry.c
	${AVALON_SRC_PATH}/forio.c
	${AVALON_SRC_PATH}/depictutil.c
	${AVALON_SRC_PATH}/denormal.c
	${AVALON_SRC_PATH}/casutils.c
	${AVALON_SRC_PATH}/ssmatch.c 
	${AVALON_SRC_PATH}/rtutils.c
	${AVALON_SRC_PATH}/smi2mol.c
	${AVALON_SRC_PATH}/didepict.c
	${AVALON_SRC_PATH}/pattern.c
	${AVALON_SRC_PATH}/canonizer.c
	${AVALON_SRC_PATH}/aacheck.c
	${AVALON_SRC_PATH}/fixcharges.c
	${AVALON_SRC_PATH}/struchk.c
	${AVALON_SRC_PATH}/reaccsio.c 
	${AVALON_SRC_PATH}/hashcode.c
  )
rdkit_library(avalon_clib ${avalon_clib_srcs})


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${AVALON_SRC_PATH})

rdkit_library(AvalonLib AvalonTools.cpp SHARED LINK_LIBRARIES avalon_clib FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral )
rdkit_test(testAvalonLib1 test1.cpp
           LINK_LIBRARIES AvalonLib avalon_clib FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral)

add_subdirectory(Wrap)

endif(RDK_BUILD_AVALON_SUPPORT)
