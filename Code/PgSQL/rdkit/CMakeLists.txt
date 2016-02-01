set(EXTENSION "rdkit")
include_directories(${PostgreSQL_INCLUDE_DIRS})
if(WIN32)
  include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32)
  if(MSVC)
    include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32_msvc)
  endif(MSVC)
endif(WIN32)
if(RDK_BUILD_AVALON_SUPPORT)
  include_directories(${CMAKE_SOURCE_DIR}/External)
  add_definitions(-DBUILD_AVALON_SUPPORT)
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  include_directories(${CMAKE_SOURCE_DIR}/External)
  add_definitions(-DBUILD_INCHI_SUPPORT)
endif(RDK_BUILD_INCHI_SUPPORT)
add_definitions(-DRDKITVER="007300")
if(MSVC)
  add_definitions(-DWIN32 -DBUILDING_MODULE -DNOMINMAX)
endif(MSVC)
link_directories(${PostgreSQL_LIBRARY_DIRS})
add_library(${EXTENSION} SHARED
            adapter.cpp bfp_op.c cache.c guc.c low_gist.c mol_op.c
            rdkit_gist.c rdkit_io.c rxn_op.c sfp_op.c)
target_link_libraries(${EXTENSION} ${PostgreSQL_LIBRARIES})
if(WIN32)
  target_link_libraries(${EXTENSION} postgres)
endif(WIN32)
set(avalonRegress "")
set(inchiRegress "")
if(RDK_BUILD_AVALON_SUPPORT)
  set(avalonRegress "avalon")
  target_link_libraries(${EXTENSION} AvalonLib avalon_clib)
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  set(inchiRegress "inchi")
  target_link_libraries(${EXTENSION} RDInchiLib Inchi)
endif(RDK_BUILD_INCHI_SUPPORT)
target_link_libraries(${EXTENSION} MolHash FMCS ChemReactions ChemTransforms
                      FileParsers SmilesParse Fingerprints Subgraphs Descriptors
                      PartialCharges SubstructMatch GraphMol DataStructs Depictor
                      RDGeometryLib RDGeneral)
if(RDK_BUILD_THREADSAFE_SSS)
  target_link_libraries(${EXTENSION} ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_THREADSAFE_SSS)

if(NOT DEFINED ${PostgreSQL_ROOT})
  if(DEFINED ENV{PostgreSQL_ROOT})
    set(PostgreSQL_ROOT "$ENV{PostgreSQL_ROOT}")
  endif(DEFINED ENV{PostgreSQL_ROOT})
endif(NOT DEFINED ${PostgreSQL_ROOT})
if(NOT DEFINED ${PostgreSQL_ROOT})
  set(PostgreSQL_ROOT "${PostgreSQL_LIBRARY_DIRS}/..")
endif(NOT DEFINED ${PostgreSQL_ROOT})
if(NOT DEFINED ${PostgreSQL_CONFIG_DIR})
  set(PostgreSQL_CONFIG_DIR "${PostgreSQL_ROOT}/bin")
endif(NOT DEFINED ${PostgreSQL_CONFIG_DIR})
macro (run_pg_config arg var)
  execute_process(COMMAND ${PostgreSQL_CONFIG_DIR}/pg_config ${arg}
                  RESULT_VARIABLE pgsql_config_result
                  OUTPUT_VARIABLE ${var}
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT ${pgsql_config_result} EQUAL 0 OR NOT ${var})
    message(FATAL_ERROR "${PostgreSQL_CONFIG_DIR}/pg_config ${arg} failed")
  endif()
endmacro ()
run_pg_config (--bindir PG_BINDIR)
run_pg_config (--sharedir PG_SHAREDIR)
run_pg_config (--pkglibdir PG_PKGLIBDIR)
set(PG_PKGLIBDIR "${PG_PKGLIBDIR}/")
set(PG_EXTENSIONDIR "${PG_SHAREDIR}/extension")
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
set(testPgSQLName "${CMAKE_CURRENT_BINARY_DIR}/pgsql_regress")
set(installPgSQLName "${CMAKE_CURRENT_BINARY_DIR}/pgsql_install")
set(PG_CURRENT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/")
set(PG_CURRENT_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/")
if(WIN32 AND NOT UNIX)
  set(testPgSQLCommand "")
  set(testPgSQLName "${testPgSQLName}.bat")
  set(installPgSQLName "${installPgSQLName}.bat")
  set(installCopyCommand "copy /Y")
  set(PG_RDKIT_LIB_SRC "${EXTENSION}.dll")
  set(PG_RDKIT_LIB_DEST "${EXTENSION}.dll")
  set(PGREGRESS_BINDIR_SWITCH "bindir")
  string(REGEX REPLACE "/" "\\\\" PG_CURRENT_SOURCE_DIR ${PG_CURRENT_SOURCE_DIR})
  string(REGEX REPLACE "/" "\\\\" PG_CURRENT_BINARY_DIR ${PG_CURRENT_BINARY_DIR})
  string(REGEX REPLACE "/" "\\\\" PG_PKGLIBDIR ${PG_PKGLIBDIR})
else(WIN32 AND NOT UNIX)
  set(testPgSQLCommand "sh")
  set(testPgSQLName "${testPgSQLName}.sh")
  set(installPgSQLName "${installPgSQLName}.sh")
  set(installCopyCommand "cp")
  set(PG_RDKIT_LIB_SRC "lib${EXTENSION}.so")
  set(PG_RDKIT_LIB_DEST "${EXTENSION}.so")
  set(PGREGRESS_BINDIR_SWITCH "psqldir")
endif(WIN32 AND NOT UNIX)
set(testPgSQLBody "cd \"${PG_CURRENT_SOURCE_DIR}\"\n"
    "\"${PGREGRESS_BINARY}\" --inputdir=sql "
    "--${PGREGRESS_BINDIR_SWITCH}=\"${PG_BINDIR}\" rdkit-91 "
    "props btree molgist bfpgist-91 sfpgist slfpgist fps reaction "
    "${avalonRegress} ${inchiRegress}\n")
file(STRINGS ${PG_CURRENT_SOURCE_DIR}${EXTENSION}.control
    PG_EXTVERSION LIMIT_COUNT 1 REGEX default_version)
string(REPLACE " " ";" tokList ${PG_EXTVERSION})
list(GET tokList -1 PG_EXTVERSION)
string(REGEX REPLACE "'" "" PG_EXTVERSION ${PG_EXTVERSION})
if(MSVC)
  set(PG_RDKIT_LIB_SRC "${CMAKE_BUILD_TYPE}\\${PG_RDKIT_LIB_SRC}")
endif(MSVC)
set(installPgSQLBody
    "${installCopyCommand} \"${PG_CURRENT_SOURCE_DIR}${EXTENSION}.sql91.in\" "
    "\"${PG_CURRENT_SOURCE_DIR}${EXTENSION}--${PG_EXTVERSION}.sql\"\n"
    "${installCopyCommand} \"${PG_CURRENT_SOURCE_DIR}${EXTENSION}--${PG_EXTVERSION}.sql\" \"${PG_EXTENSIONDIR}\"\n"
    "${installCopyCommand} \"${PG_CURRENT_SOURCE_DIR}${EXTENSION}.control\" \"${PG_EXTENSIONDIR}\"\n"
    "${installCopyCommand} \"${PG_CURRENT_BINARY_DIR}${PG_RDKIT_LIB_SRC}\" "
    "\"${PG_PKGLIBDIR}${PG_RDKIT_LIB_DEST}\"\n")
file(WRITE ${testPgSQLName} ${testPgSQLBody})
file(WRITE ${installPgSQLName} ${installPgSQLBody})
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
