cmake_minimum_required(VERSION 3.2)
project(coordgen)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options
option(COORDGEN_RIGOROUS_BUILD "Abort the build if the compiler issues \
       any warnings" ON )
option(COORDGEN_BUILD_TESTS "Whether test executables should be built" ON)
option(COORDGEN_BUILD_EXAMPLE "Whether to build the sample executable" ON)
option(COORDGEN_USE_MAEPARSER "Whether to allow loading of run-time templates" OFF)
option(COORDGEN_BUILD_SHARED_LIBS "Build coordgen as a shared library \
       (turn off for a static one)" ON)
option(DEBUG_COORDS "Whether to print minimization data and atom mapping to file" OFF)

# Use the maeparser_DIR variable to tell CMake where to search for the
# maeparser library.

if (COORDGEN_USE_MAEPARSER OR COORDGEN_BUILD_TESTS)
    set(USE_MAEPARSER ON)
else()
    set(USE_MAEPARSER OFF)
endif()

if (USE_MAEPARSER)
    set(MAEPARSER_VERSION "master" CACHE STRING "maeparser tag to build if \
        a compiled library is not found")
endif()

if(MSVC)
    # C4251 disables warnings for export STL containers as arguments
    # (returning a vector of things)
    add_definitions(/wd4251 /wd4275 /wd4996 /D_SCL_SECURE_NO_WARNINGS
                    /D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DBOOST_ALL_NO_LIB)
endif(MSVC)

if(COORDGEN_RIGOROUS_BUILD)
if(MSVC)
    add_definitions(/WX)
else(MSVC)
    add_definitions(-Wall -Wextra -Werror)
endif(MSVC)
endif(COORDGEN_RIGOROUS_BUILD)

# Source files & headers
file(GLOB SOURCES "*.cpp")

# Build Targets & Configuration -- coordgen library
if(COORDGEN_BUILD_SHARED_LIBS)
    add_library(coordgen SHARED ${SOURCES})
    target_compile_definitions(coordgen PRIVATE "IN_COORDGEN")
    set_property(TARGET coordgen PROPERTY CXX_VISIBILITY_PRESET "hidden")
else(COORDGEN_BUILD_SHARED_LIBS)
    add_library(coordgen STATIC ${SOURCES})
    if (USE_MAEPARSER)
        target_compile_definitions(coordgen PRIVATE "STATIC_MAEPARSER")
    endif(USE_MAEPARSER)
    target_compile_definitions(coordgen PRIVATE "STATIC_COORDGEN")
endif(COORDGEN_BUILD_SHARED_LIBS)

if (DEBUG_COORDS)
   target_compile_definitions(coordgen PRIVATE "DEBUG_MINIMIZATION_COORDINATES")
endif()

# Dependencies
if (USE_MAEPARSER)
    if(TARGET maeparser)
        message(STATUS "Using externally defined maeparser target to "
        "build coordgen")
    else()
        include(CoordgenUtils)
        set(MAEPARSER_BUILD_SHARED_LIBS ${COORDGEN_BUILD_SHARED_LIBS} CACHE BOOL "Library type for maeparser")
        find_or_clone_maeparser()
    endif()

    include_directories(${maeparser_INCLUDE_DIRS})
    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
endif (USE_MAEPARSER)

if (COORDGEN_USE_MAEPARSER)
    target_link_libraries(coordgen ${maeparser_LIBRARIES})
    target_compile_definitions(coordgen PRIVATE "USE_MAEPARSER")
endif()

set_target_properties(coordgen
    PROPERTIES
        VERSION 3.0.2
        SOVERSION 3
)

if(MSVC)
    set(CMAKE_INSTALL_LIBDIR lib)
    set(CMAKE_INSTALL_BINDIR bin)
else(MSVC)
    include(GNUInstallDirs)
endif(MSVC)

# Install configuration
install(TARGETS coordgen
    EXPORT coordgen-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES
    CoordgenConfig.hpp
    CoordgenFragmentBuilder.h
    CoordgenFragmenter.h
    CoordgenMacrocycleBuilder.h
    CoordgenMinimizer.h
    sketcherMinimizerAtom.h
    sketcherMinimizerBendInteraction.h
    sketcherMinimizerBond.h
    sketcherMinimizerClashInteraction.h
    sketcherMinimizerEZConstrainInteraction.h
    sketcherMinimizerFragment.h
    sketcherMinimizer.h
    sketcherMinimizerInteraction.h
    sketcherMinimizerMarchingSquares.h
    sketcherMinimizerMaths.h
    sketcherMinimizerMolecule.h
    sketcherMinimizerResidue.h
    sketcherMinimizerResidueInteraction.h
    sketcherMinimizerRing.h
    sketcherMinimizerStretchInteraction.h
    DESTINATION include/coordgen)

install(EXPORT coordgen-targets
    FILE ${PROJECT_NAME}-config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)

# Example
if(COORDGEN_BUILD_EXAMPLE)
    add_subdirectory(example_dir)
endif(COORDGEN_BUILD_EXAMPLE)

# Tests
if(COORDGEN_BUILD_TESTS)
    set(MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --time-stamp=yes \
        --num-callers=20 --gen-suppressions=all --leak-check=full \
        --show-reachable=no --trace-children=yes --error-exitcode=29")
    include(CTest)
    add_subdirectory(test)
endif(COORDGEN_BUILD_TESTS)
