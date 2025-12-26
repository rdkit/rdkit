remove_definitions(-DRDKIT_RDGEOMETRYLIB_BUILD)

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
# if("${PYTHON_LDSHARED}" STREQUAL "")
# else()
#    set_target_properties(rd PROPERTIES LINK_FLAGS ${PYTHON_LDSHARED})
#endif()
set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")

nanobind_add_module(rdGeometry NB_SHARED rdGeometry.cpp)
target_link_libraries(rdGeometry PUBLIC RDGeometryLib rdkit_base)

nanobind_add_module(rdGeometryTest NB_SHARED rdGeometryTest.cpp)
target_link_libraries(rdGeometryTest PUBLIC RDGeometryLib rdkit_base)


INSTALL(TARGETS rdGeometry rdGeometryTest
        LIBRARY DESTINATION ${RDKit_PythonDir}/Geometry COMPONENT python)

set_target_properties(nanobind
  PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY  ${RDK_LIBRARY_OUTPUT_DIRECTORY} )



#rdkit_python_extension(rdGeometry
#                       Point.cpp UniformGrid3D.cpp UniformRealValueGrid3D.cpp rdGeometry.cpp
#                       DEST Geometry
#                       LINK_LIBRARIES RDGeometryLib)

#add_pytest(testPyGeometry ${CMAKE_CURRENT_SOURCE_DIR}/testGeometry.py)

