project (jRDKFuncs Java)

find_package(Java REQUIRED)
find_package(JNI REQUIRED)
INCLUDE_DIRECTORIES(${JNI_INCLUDE_DIRS})

SET_SOURCE_FILES_PROPERTIES(RDKFuncs.i PROPERTIES CPLUSPLUS ON )

SET(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/org/RDKit )
SET(CMAKE_SWIG_FLAGS -package "org.RDKit")
SET(SWIG_MODULE_RDKFuncs_EXTRA_DEPS ../RDKFuncs.i)
SWIG_ADD_MODULE(RDKFuncs "java" RDKFuncs.i ../RDKFuncs.cpp)

SWIG_LINK_LIBRARIES(RDKFuncs ${RDKit_Wrapper_Libs} )

# code adapted from the wrapper code for
# GDCM: http://gdcm.svn.sf.net/viewvc/gdcm/trunk/Wrapping/Java/CMakeLists.txt?view=markup
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar
## 1. run this custom command only after swig has been run.
  COMMAND ${JAVA_COMPILE} ARGS "*.java" "org/RDKit/*.java" 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
## 2. now that the *.class have been generated construct the jar file.
  COMMAND ${JAVA_ARCHIVE} ARGS cf 
    ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar
    org/RDKit/*.class 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS "${swig_generated_file_fullname}"
  COMMENT "javac *.java; jar cf -> .jar"
)
ADD_CUSTOM_TARGET(RDKitJava ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar
  COMMENT "building jar"
)


INSTALL(TARGETS RDKFuncs 
              DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
