project (jGraphMolJava Java)

include_directories( ${RDKit_ExternalDir} )

find_package(Java REQUIRED)
find_package(JNI REQUIRED)

# Make sure we can access the junit.jar -- needed for compiling/running tests
find_file(JUNIT_JAR junit.jar ${RDKit_JavaLibDir})
if (NOT JUNIT_JAR)
	MESSAGE ("junit.jar is not found in ${RDKit_JavaLibDir}.  Please add it to this directory and rerun cmake.")
	MESSAGE(FATAL_ERROR "Cannot find required JUnit library (junit.jar)")
endif (NOT JUNIT_JAR)

# Likewise for javadoc
find_program(JAVADOC_EXE javadoc PATH _JAVA_PATHS _JAVA_HINTS /usr/bin /usr/local/bin)
if (NOT JAVADOC_EXE)
	MESSAGE ("javadoc (executable) is not found. Please add it to PATH and rerun cmake.")
	MESSAGE(FATAL_ERROR "Cannot find required javadoc executable (javadoc)")
endif (NOT JAVADOC_EXE)

INCLUDE_DIRECTORIES(${JNI_INCLUDE_DIRS})

SET_SOURCE_FILES_PROPERTIES(GraphMolJava.i PROPERTIES CPLUSPLUS ON )

FILE(TO_NATIVE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../*.java COPY_SOURCE)
FILE(TO_NATIVE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/org/RDKit COPY_DEST)

# Setup a few variables for environment-specific things
if(WIN32)
  ADD_DEFINITIONS("/W3 /wd4716 /bigobj")
  SET(PATH_SEP ";")
  SET(COPY_CMD xcopy ${COPY_SOURCE} ${COPY_DEST} /Y /I)
else()
  SET(PATH_SEP ":")
  SET(COPY_CMD cp -p ${COPY_SOURCE} ${COPY_DEST})
endif()

# Coax SWIG into playing nicely with Apple environments
if(APPLE)
  SET(CMAKE_SIZEOF_VOID_P 4)
endif(APPLE)

if(CMAKE_SIZEOF_VOID_P MATCHES 4)
  SET(CMAKE_SWIG_FLAGS -package "org.RDKit")
else()
  SET(CMAKE_SWIG_FLAGS -package "org.RDKit" "-DSWIGWORDSIZE64")
endif()
SET(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/org/RDKit )

if(RDK_BUILD_INCHI_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DBUILD_INCHI_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()
if(RDK_BUILD_AVALON_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DBUILD_AVALON_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()


# enable this line to build the ErrorGenerator class for testing handling of C++ errors in the JNI layer
#SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-D INCLUDE_ERROR_GENERATOR" )

FILE(GLOB SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../*.i")

# we added all source files, now remove the ones that we're not supporting in this build:
if(NOT RDK_BUILD_AVALON_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../AvalonLib.i")
endif()

if(NOT RDK_BUILD_INCHI_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../Inchi.i")
endif()

SET(SWIG_MODULE_GraphMolWrap_EXTRA_DEPS ${SWIG_SRC_FILES} )

SWIG_ADD_MODULE(GraphMolWrap "java" GraphMolJava.i )

# it doesnt seem like the threading libs should need to be here, but
# as of Oct 2012 using boost 1.51 under at least ubuntu 12.04 we get a
# link error if they aren't there.
SWIG_LINK_LIBRARIES(GraphMolWrap ${RDKit_Wrapper_Libs}
      ${RDKit_THREAD_LIBS})

# code adapted from the wrapper code for
# GDCM: http://gdcm.svn.sf.net/viewvc/gdcm/trunk/Wrapping/Java/CMakeLists.txt?view=markup
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar

  COMMAND ${CMAKE_COMMAND} -E make_directory build
  COMMAND ${CMAKE_COMMAND} -E make_directory build-test

  ## Add in our own Java sources
  #COMMAND ${COPY_CMD}
  ## 1. run this custom command only after swig has been run.
  COMMAND ${JAVA_COMPILE} -d build "src/org/RDKit/*.java"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  ## 2. now that the *.class have been generated construct the jar file.
  COMMAND ${JAVA_ARCHIVE} cf
    ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar -C ${CMAKE_CURRENT_SOURCE_DIR}/build
    org/RDKit
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  #DEPENDS GraphMolWrap
  DEPENDS "${swig_generated_file_fullname}"
  COMMENT "javac *.java; jar cf -> .jar"

  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKitDoc.jar
  ## Run Javadoc against newly created .java files to create docs
  COMMAND ${JAVADOC_EXE} -tag notes -tag example -d ${CMAKE_CURRENT_SOURCE_DIR}/doc -sourcepath ${CMAKE_CURRENT_SOURCE_DIR}/src org.RDKit
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  ## Put the doc files into their own separate archive.
  COMMAND ${JAVA_ARCHIVE} cf 
    ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKitDoc.jar
    doc
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  #DEPENDS GraphMolWrap
  DEPENDS "${swig_generated_file_fullname}"
  COMMENT "jar cf org.RDKitDoc.jar -d doc"
)
ADD_CUSTOM_TARGET(GraphMolWrapJar ALL
  DEPENDS GraphMolWrap ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar
  COMMENT "building jar"
)

## Tests -- note that building/rebuilding them keys on the file WrapperTests.class

SET(CMAKE_JAVA_TEST_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/build-test )
FILE(GLOB JAVA_TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src-test/org/RDKit/*.java")

# we added all source files, now remove the ones that we're not supporting in this build:
if(NOT RDK_BUILD_AVALON_SUPPORT)
LIST(REMOVE_ITEM JAVA_TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src-test/org/RDKit/AvalonLibTests.java")
endif()

if(NOT RDK_BUILD_INCHI_SUPPORT)
LIST(REMOVE_ITEM JAVA_TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src-test/org/RDKit/InchiTests.java")
endif()

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_JAVA_TEST_OUTDIR}/org/RDKit/WrapperTests.class
  COMMAND ${JAVA_COMPILE} -d ${CMAKE_JAVA_TEST_OUTDIR} -cp "${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar${PATH_SEP}${JUNIT_JAR}" ${JAVA_TEST_FILES}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar ${JAVA_TEST_FILES}
)
ADD_CUSTOM_TARGET(BuildJavaWrapperTests ALL
  DEPENDS ${CMAKE_JAVA_TEST_OUTDIR}/org/RDKit/WrapperTests.class GraphMolWrap 
  COMMENT "building test classes"
)
ADD_CUSTOM_TARGET(RunJavaWrapperTests 
  DEPENDS BuildJavaWrapperTests
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runUnitTests.sh
)

ADD_TEST(JavaAromaticTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.AromaticTests)
ADD_TEST(JavaAtomPairsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.AtomPairsTests)
ADD_TEST(JavaBasicMoleculeTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.BasicMoleculeTests)
ADD_TEST(JavaBasicMolecule2Tests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.BasicMolecule2Tests)
ADD_TEST(JavaChemAtomTests
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ChemAtomTests)
ADD_TEST(JavaChemBondTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ChemBondTests)
ADD_TEST(JavaChemReactionTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.ChemReactionTests)
ADD_TEST(JavaChemSmartsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ChemSmartsTests)
ADD_TEST(JavaChemTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ChemTests)
ADD_TEST(JavaChemv2Tests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.Chemv2Tests)
ADD_TEST(JavaConformerTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ConformerTests)
ADD_TEST(JavaDescriptorTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.DescriptorTests)
ADD_TEST(JavaDistanceGeometryTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.DistanceGeometryTests)
ADD_TEST(JavaErrorHandlingTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ErrorHandlingTests)
ADD_TEST(JavaFingerprintsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.FingerprintsTests)
ADD_TEST(JavaForceFieldsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.ForceFieldsTests)
ADD_TEST(JavaHManipulationsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.HManipulationsTests)
ADD_TEST(JavaLipinskiTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.LipinskiTests)
ADD_TEST(JavaPicklingTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.PicklingTests)
ADD_TEST(JavaSmilesCreationTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.SmilesCreationTests)
ADD_TEST(JavaSmilesDetailsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.SmilesDetailsTests)
ADD_TEST(JavaSmilesTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.SmilesTests)
ADD_TEST(JavaSuppliersTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.SuppliersTests)
ADD_TEST(JavaWrapperTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.WrapperTests)

if(RDK_BUILD_AVALON_SUPPORT)
   ADD_TEST(JavaAvalonTests 
      java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
  	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
	org.RDKit.AvalonLibTests)
endif(RDK_BUILD_AVALON_SUPPORT)       

if (RDK_BUILD_INCHI_SUPPORT)
   ADD_TEST(JavaInchiTests 
       java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR}
       	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
		org.RDKit.InchiTests)
endif (RDK_BUILD_INCHI_SUPPORT)

#ADD_TEST(JavaMemoryTests 
#    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
#	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar" 
#	org.RDKit.MemoryTests)

ADD_TEST(JavaChemTransformsTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.ChemTransformsTests)
ADD_TEST(JavaFMCSTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.FMCSTests)
ADD_TEST(JavaPDBTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.PDBTests)

ADD_TEST(JavaAlignTests 
    java -Djava.library.path=${CMAKE_CURRENT_SOURCE_DIR} 
	-cp "${JUNIT_JAR}${PATH_SEP}${CMAKE_JAVA_TEST_OUTDIR}${PATH_SEP}${CMAKE_CURRENT_SOURCE_DIR}/org.RDKit.jar"
	org.RDKit.AlignTests)

INSTALL(TARGETS GraphMolWrap 
              DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
