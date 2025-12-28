remove_definitions(-DRDKIT_RDGEOMETRYLIB_BUILD)


# -------------
#  This block is stuff that should land in the main CMakeLists.txt file

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

find_package(Python 3.9 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)
set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")
#----------------

# Should add an rdkit function for this and the INSTALL.
nanobind_add_module(rdGeometry NB_SHARED rdGeometry.cpp)
target_link_libraries(rdGeometry PUBLIC RDGeometryLib rdkit_base)

nanobind_add_module(rdGeometryTest NB_SHARED rdGeometryTest.cpp)
target_link_libraries(rdGeometryTest PUBLIC RDGeometryLib rdkit_base)

INSTALL(TARGETS rdGeometry rdGeometryTest
        LIBRARY DESTINATION ${RDKit_PythonDir}/Geometry COMPONENT python)

#-------
#  these also should go somewhere else, but they
#    need to appear afer the first call to nanobind_add_module
set_target_properties(nanobind
  PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY  ${RDK_LIBRARY_OUTPUT_DIRECTORY} )
#INSTALL(TARGETS nanobind LIBRARY DESTINATION ${RDK_LIBRARY_OUTPUT_DIRECTORY})

INSTALL(TARGETS nanobind 
        DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
        COMPONENT ${sharedLibComponent})


#-------

#add_pytest(testPyGeometry ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py)

