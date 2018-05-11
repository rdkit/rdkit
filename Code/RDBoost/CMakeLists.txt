if(WIN32 OR "${Py_ENABLE_SHARED}" STREQUAL "1")

  add_definitions(-DRDKIT_RDBOOST_BUILD)
  rdkit_library(RDBoost Wrap.cpp
                LINK_LIBRARIES RDGeneral ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
else()

    add_definitions(-DRDKIT_RDBOOST_BUILD)
    rdkit_library(RDBoost Wrap.cpp
                LINK_LIBRARIES RDGeneral ${Boost_LIBRARIES})
    if("${PYTHON_LDSHARED}" STREQUAL "")
    else()
       # fixes an apple issue
       STRING(REGEX REPLACE "-bundle" "" LDSHARED ${PYTHON_LDSHARED})
       set_target_properties(RDBoost PROPERTIES LINK_FLAGS ${LDSHARED})
    endif()
endif()

rdkit_headers(Wrap.h PySequenceHolder.h
              list_indexing_suite.hpp
              python.h
              boost_numpy.h
              python_streambuf.h
              DEST RDBoost)

add_subdirectory(Wrap)

if(RDK_BUILD_PYTHON_WRAPPERS)
  # this is a particularly nasty little bit of hackery.
  # cmake (as of 2.8) ignores the ADDITIONAL_MAKE_CLEAN_FILES argument
  # unless it's in a directory that actualy has a target. So in order to
  # cleanup files created by the python tests we have to put something
  # here with the ugly relative path
  IF(COMMAND cmake_policy)
     CMAKE_POLICY(SET CMP0009 NEW)
  ENDIF(COMMAND cmake_policy)
  FILE(GLOB_RECURSE RDK_PYC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../../rdkit/*.pyc")
  FILE(GLOB_RECURSE RDK_PROJECTSQLT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../../Projects/*.sqlt")
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
        "${RDK_PYC_FILES};${RDK_DATASQLT_FILES};${RDK_PROJECTSQLT_FILES}")
endif(RDK_BUILD_PYTHON_WRAPPERS)
