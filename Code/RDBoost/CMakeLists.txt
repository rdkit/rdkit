rdkit_library(RDBoost Wrap.cpp
              LINK_LIBRARIES ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

rdkit_headers(Exceptions.h Wrap.h PySequenceHolder.h
              list_indexing_suite.hpp
              python_streambuf.h
              DEST RDBoost)

add_subdirectory(Wrap)

if(RDK_BUILD_PYTHON_WRAPPERS)
  # this is a particularly nasty little bit of hackery.
  # cmake (as of 2.8) ignores the ADDITIONAL_MAKE_CLEAN_FILES argument
  # unless it's in a directory that actualy has a target. So in order to 
  # cleanup files created by the python tests we have to put something
  # here with the ugly relative path
  FILE(GLOB_RECURSE RDK_PYC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../../rdkit/*.pyc")
  FILE(GLOB_RECURSE RDK_PROJECTSQLT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../../Projects/*.sqlt")
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
        "${RDK_PYC_FILES};${RDK_DATASQLT_FILES};${RDK_PROJECTSQLT_FILES}")
endif(RDK_BUILD_PYTHON_WRAPPERS)
