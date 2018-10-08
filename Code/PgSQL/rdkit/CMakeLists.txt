cmake_minimum_required(VERSION 3.1)

if(APPLE)
  set (EXTENSION_PREFIX "")
  set (EXTENSION_SUFFIX ".so")
else(APPLE)
  set (EXTENSION_PREFIX "lib")
  set (EXTENSION_SUFFIX "")
endif(APPLE)

find_package(PostgreSQL REQUIRED)

if(RDK_OPTIMIZE_NATIVE)
  add_definitions(-DUSE_BUILTIN_POPCOUNT)
endif(RDK_OPTIMIZE_NATIVE)

set(EXTENSION rdkit)
if(NOT DEFINED Boost_INCLUDE_DIRS)
  find_package(Boost 1.56.0 REQUIRED)
endif(NOT DEFINED Boost_INCLUDE_DIRS)
include_directories(${Boost_INCLUDE_DIRS})
message("postgres: ${PostgreSQL_INCLUDE_DIRS}")
include_directories(${PostgreSQL_INCLUDE_DIRS})
if(WIN32)
  include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32)
  if(MSVC)
    include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32_msvc)
  endif(MSVC)
endif(WIN32)
if(RDK_BUILD_AVALON_SUPPORT)
  include_directories(${CMAKE_SOURCE_DIR}/External)
  add_definitions(-DRDK_BUILD_AVALON_SUPPORT)
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  include_directories(${CMAKE_SOURCE_DIR}/External)
  add_definitions(-DRDK_BUILD_INCHI_SUPPORT)
endif(RDK_BUILD_INCHI_SUPPORT)
add_definitions(-DRDKITVER="007400")
if(NOT DEFINED PostgreSQL_ROOT)
  if(DEFINED ENV{PostgreSQL_ROOT})
    set(PostgreSQL_ROOT "$ENV{PostgreSQL_ROOT}")
  endif(DEFINED ENV{PostgreSQL_ROOT})
endif(NOT DEFINED PostgreSQL_ROOT)
if(NOT DEFINED PostgreSQL_ROOT)
  set(PostgreSQL_ROOT "${PostgreSQL_LIBRARY_DIRS}/..")
endif(NOT DEFINED PostgreSQL_ROOT)
if(NOT DEFINED PostgreSQL_CONFIG_DIR)
  set(PostgreSQL_CONFIG_DIR "${PostgreSQL_ROOT}/bin")
endif(NOT DEFINED PostgreSQL_CONFIG_DIR)
if(NOT DEFINED PostgreSQL_CONFIG)
  set(PostgreSQL_CONFIG "${PostgreSQL_CONFIG_DIR}/pg_config")
endif(NOT DEFINED PostgreSQL_CONFIG)
macro (run_pg_config arg var)
  execute_process(COMMAND ${PostgreSQL_CONFIG} ${arg}
                  RESULT_VARIABLE pgsql_config_result
                  OUTPUT_VARIABLE ${var}
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT ${pgsql_config_result} EQUAL 0 OR NOT ${var})
    message(FATAL_ERROR "${PostgreSQL_CONFIG} ${arg} failed")
  endif()
endmacro ()
run_pg_config (--bindir PG_BINDIR)
run_pg_config (--sharedir PG_SHAREDIR)
run_pg_config (--pkglibdir PG_PKGLIBDIR)
set(PG_PKGLIBDIR "${PG_PKGLIBDIR}/")
set(PG_EXTENSIONDIR "${PG_SHAREDIR}/extension")
if(MSVC)
  add_definitions(-DWIN32 -DBUILDING_MODULE -DNOMINMAX)
  if(NOT (MSVC_VERSION LESS 1700))
    add_definitions(-DHAVE_RINT)
  endif()
endif(MSVC)
if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "-bundle -multiply_defined suppress"
              "-Wl,-dead_strip_dylibs -bundle_loader ${PG_BINDIR}/postgres")
  add_executable("${EXTENSION}${EXTENSION_SUFFIX}"
              adapter.cpp bfp_op.c cache.c guc.c low_gist.c mol_op.c
              rdkit_gist.c bfp_gist.c bfp_gin.c bitstring.c rdkit_io.c rxn_op.c sfp_op.c)
else(APPLE)
  link_directories(${PostgreSQL_LIBRARY_DIRS})
  link_directories(${RDKit_LibDir})
  add_library(${EXTENSION} SHARED
              adapter.cpp bfp_op.c cache.c guc.c low_gist.c mol_op.c
              rdkit_gist.c bfp_gist.c bfp_gin.c bitstring.c rdkit_io.c rxn_op.c sfp_op.c)
  target_link_libraries(${EXTENSION} ${PostgreSQL_LIBRARIES})
  if(WIN32)
    target_link_libraries(${EXTENSION} postgres)
  endif(WIN32)
endif(APPLE)
set(avalonRegress "")
set(inchiRegress "")
set(pgRDKitLibSuffix "")
if(RDK_PGSQL_STATIC AND ((NOT MSVC) OR (MSVC AND RDK_INSTALL_DLLS_MSVC)))
  set(pgRDKitLibSuffix "_static")
endif()
set(pgRDKitLibList "")
if(RDK_BUILD_AVALON_SUPPORT)
  set(avalonRegress "avalon")
  set(pgRDKitLibList "${pgRDKitLibList}AvalonLib;avalon_clib;")
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  if(NOT DEFINED INCHI_LIBRARIES)
    set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/")
    if(NOT DEFINED RDKit_RELEASENAME)
      include(RDKitUtils)
    endif(NOT DEFINED RDKit_RELEASENAME)
    find_package(Inchi)
  endif(NOT DEFINED INCHI_LIBRARIES)
  set(inchiRegress "inchi")
  set(pgRDKitLibList "${pgRDKitLibList}RDInchiLib;${INCHI_LIBRARIES};")
endif(RDK_BUILD_INCHI_SUPPORT)
set(pgRDKitLibList "${pgRDKitLibList}"
    "MolDraw2D;MolTransforms;MolHash;FMCS;ChemReactions;ChemTransforms;"
    "FileParsers;SmilesParse;Fingerprints;Subgraphs;Descriptors;"
    "PartialCharges;SubstructMatch;GraphMol;EigenSolvers;DataStructs;Depictor;"
    "RDGeometryLib;RDGeneral")
foreach(pgRDKitLib ${pgRDKitLibList})
  set(pgRDKitLibs "${pgRDKitLibs}${pgRDKitLib}${pgRDKitLibSuffix};")
endforeach()
target_link_libraries(${EXTENSION}${EXTENSION_SUFFIX} ${pgRDKitLibs})
if(RDK_BUILD_THREADSAFE_SSS)
  target_link_libraries(${EXTENSION}${EXTENSION_SUFFIX} ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_THREADSAFE_SSS)

if(WIN32 AND NOT UNIX)
  set(PGREGRESS_BINARY "${PG_BINDIR}\\pg_regress.exe")
  string(REGEX REPLACE "/" "\\\\" PG_EXTENSIONDIR ${PG_EXTENSIONDIR})
  string(REGEX REPLACE "/" "\\\\" PGREGRESS_BINARY ${PGREGRESS_BINARY})
else(WIN32 AND NOT UNIX)
  run_pg_config (--pgxs PG_MAKEFILE)
  get_filename_component(PG_MAKEFILESDIR ${PG_MAKEFILE} PATH)
  set(PGREGRESS_BINARY "${PG_MAKEFILESDIR}/../test/regress/pg_regress")
endif(WIN32 AND NOT UNIX)
if(NOT EXISTS ${PGREGRESS_BINARY})
  message(FATAL_ERROR "${PGREGRESS_BINARY} does not exist")
endif()

macro (run_pgregress arg var)
  execute_process(COMMAND ${PGREGRESS_BINARY} ${arg}
                  OUTPUT_VARIABLE ${var}
                  ERROR_VARIABLE ${var}
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
endmacro ()
set(PGREGRESS_BINDIR_SWITCH "")
foreach(opt "bindir" "psqldir")
  run_pgregress("--${opt}" err)
  if (${err} MATCHES "option requires an argument")
    set(PGREGRESS_BINDIR_SWITCH ${opt})
    break()
  endif()
endforeach()
if (NOT ${PGREGRESS_BINDIR_SWITCH} STREQUAL "")
  set(PGREGRESS_BINDIR_SWITCH "--${PGREGRESS_BINDIR_SWITCH}=\"${PG_BINDIR}\"")
endif()

set(testPgSQLName "${CMAKE_CURRENT_BINARY_DIR}/pgsql_regress")
set(installPgSQLName "${CMAKE_CURRENT_BINARY_DIR}/pgsql_install")
set(PG_CURRENT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/")
set(PG_CURRENT_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/")
if(WIN32 AND NOT UNIX)
  set(testPgSQLCommand "")
  set(testPgSQLName "${testPgSQLName}.bat")
  set(installPgSQLName "${installPgSQLName}.bat")
  set(installCopyCommand "copy /Y")
  set(installPgSQLBody "")
  set(PG_RDKIT_LIB_SRC "${EXTENSION}.dll")
  set(PG_RDKIT_LIB_DEST "${EXTENSION}.dll")
  string(REGEX REPLACE "/" "\\\\" PG_CURRENT_SOURCE_DIR ${PG_CURRENT_SOURCE_DIR})
  string(REGEX REPLACE "/" "\\\\" PG_CURRENT_BINARY_DIR ${PG_CURRENT_BINARY_DIR})
  string(REGEX REPLACE "/" "\\\\" PG_PKGLIBDIR ${PG_PKGLIBDIR})
else(WIN32 AND NOT UNIX)
  set(testPgSQLCommand "sh")
  set(testPgSQLName "${testPgSQLName}.sh")
  set(installPgSQLName "${installPgSQLName}.sh")
  set(installCopyCommand "cp")
  set(installPgSQLBody "set -x\n")
  set(PG_RDKIT_LIB_SRC "${EXTENSION_PREFIX}${EXTENSION}.so")
  set(PG_RDKIT_LIB_DEST "${EXTENSION}.so")
endif(WIN32 AND NOT UNIX)
set(testPgSQLBody "cd \"${PG_CURRENT_BINARY_DIR}\"\n"
    "\"${PGREGRESS_BINARY}\" --inputdir=sql "
    "${PGREGRESS_BINDIR_SWITCH} rdkit-91 "
    "props btree molgist bfpgist-91 sfpgist slfpgist fps reaction "
    "${avalonRegress} ${inchiRegress}\n")
file(STRINGS ${PG_CURRENT_SOURCE_DIR}${EXTENSION}.control
    PG_EXTVERSION LIMIT_COUNT 1 REGEX default_version)
string(REPLACE " " ";" tokList ${PG_EXTVERSION})
list(GET tokList -1 PG_EXTVERSION)
string(REGEX REPLACE "'" "" PG_EXTVERSION ${PG_EXTVERSION})
if(MSVC)
  if(${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
    set(PG_RDKIT_LIB_SRC "${PG_RDKIT_LIB_SRC}")
  else(${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
    set(PG_RDKIT_LIB_SRC "${CMAKE_BUILD_TYPE}\\${PG_RDKIT_LIB_SRC}")
  endif(${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")
endif(MSVC)
if(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.5)
    set(RDKIT_PG_BFP_GIST_FETCH "")
else(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.5)
    set(RDKIT_PG_BFP_GIST_FETCH "FUNCTION    9   gbfp_fetch (internal),")
endif(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.5)

if(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.3)
  set(RDKIT_RDK_GIN_BFP_TRICONSISTENT "")
else(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.3)
  set(RDKIT_GIN_BFP_TRICONSISTENT "FUNCTION    6   gin_bfp_triconsistent(internal, int2, bfp, int4, internal, internal, internal),")
endif(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.3)

if(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.6)
    set(RDKIT_PARALLEL_SAFE "")
else(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.6)
    set(RDKIT_PARALLEL_SAFE "PARALLEL SAFE")
endif(${PostgreSQL_VERSION_STRING} VERSION_LESS 9.6)


configure_file("${PG_CURRENT_SOURCE_DIR}${EXTENSION}.sql.in"
               "${PG_CURRENT_BINARY_DIR}${EXTENSION}--${PG_EXTVERSION}.sql")
set(installPgSQLBody "${installPgSQLBody}"
    "${installCopyCommand} \"${PG_CURRENT_BINARY_DIR}${EXTENSION}--${PG_EXTVERSION}.sql\" \"${PG_EXTENSIONDIR}\"\n"
    "${installCopyCommand} \"${PG_CURRENT_SOURCE_DIR}${EXTENSION}.control\" \"${PG_EXTENSIONDIR}\"\n"
    "${installCopyCommand} \"${PG_CURRENT_BINARY_DIR}${PG_RDKIT_LIB_SRC}\" "
    "\"${PG_PKGLIBDIR}${PG_RDKIT_LIB_DEST}\"\n")
file(WRITE ${testPgSQLName} ${testPgSQLBody})
file(WRITE ${installPgSQLName} ${installPgSQLBody})
foreach(testDir sql data expected)
  add_custom_command(TARGET ${EXTENSION}${EXTENSION_SUFFIX} PRE_BUILD
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
      "${PG_CURRENT_SOURCE_DIR}${testDir}" "${PG_CURRENT_BINARY_DIR}${testDir}")
endforeach()
message("\n"
    "=====================================================================\n"
    "After building, check\n"
    "${installPgSQLName}\n"
    "for correctness of installation paths. If everything is OK, gain \n"
    "administrator privileges, stop the PostgreSQL service, run\n"
    "${installPgSQLName}\n"
    "to install the PostgreSQL RDKit cartridge, then start again\n"
    "the PostgreSQL service\n"
    "=====================================================================\n")
add_test(testPgSQL ${testPgSQLCommand} ${testPgSQLName})
