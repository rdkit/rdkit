include_directories(${PostgreSQL_INCLUDE_DIR}
              ${PostgreSQL_INCLUDE_DIR}/server)
if(WIN32)
  include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32)
endif(WIN32)
if(MSVC)
  include_directories(${PostgreSQL_INCLUDE_DIR}/server/port/win32_msvc)
endif(MSVC)
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
link_directories(${PostgreSQL_LIBRARY_DIR})
add_library(rdkit SHARED
              adapter.cpp bfp_op.c cache.c guc.c low_gist.c mol_op.c
              rdkit_gist.c rdkit_io.c rxn_op.c sfp_op.c)
target_link_libraries(rdkit postgres)
set(avalonRegress "")
set(inchiRegress "")
if(RDK_BUILD_AVALON_SUPPORT)
  set(avalonRegress "avalon")
  target_link_libraries(rdkit AvalonLib avalon_clib)
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  set(inchiRegress "inchi")
  target_link_libraries(rdkit RDInchiLib Inchi)
endif(RDK_BUILD_INCHI_SUPPORT)
target_link_libraries(rdkit MolHash FMCS ChemReactions ChemTransforms
              FileParsers SmilesParse Fingerprints Subgraphs Descriptors
              PartialCharges SubstructMatch GraphMol DataStructs Depictor
              RDGeometryLib RDGeneral)
if(RDK_BUILD_THREADSAFE_SSS)
  target_link_libraries(rdkit ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_THREADSAFE_SSS)

# The test is a bit involved just because the much simpler "add_test"
# command which is currently commented out (see below) requires
# "-C Release" to be passed along "ctest", which is inconvenient.
# Hence, since to avoid this inconvenience we use the basic "add_test"
# syntax which does not allow to pass a WORKING_DIRECTORY parameter,
# and therefore we need to write a small bash/bat script
# to cd in the correct directory upfront. A bit clumsy, but it works
set(PGREGRESS_BINARY ${PostgreSQL_ROOT})
STRING(REGEX REPLACE "\\\\" "/" PGREGRESS_BINDIR ${PostgreSQL_ROOT})
if(WIN32 AND NOT UNIX)
  set(PGREGRESS_BINARY "${PGREGRESS_BINARY}\\bin\\pg_regress.exe")
else(WIN32 AND NOT UNIX)
  set(PGREGRESS_BINARY "${PGREGRESS_BINDIR}/bin/pg_regress")
endif(WIN32 AND NOT UNIX)
set(testPgSQLName "${CMAKE_SOURCE_DIR}/Code/PgSQL/rdkit/pg_regress")
set(PGREGRESS_BINDIR "${PGREGRESS_BINDIR}/bin")
if(WIN32 AND NOT UNIX)
  set(testPgSQLName "${testPgSQLName}.bat")
  set(testPgSQLHead "cd %~dp0\n")
else(WIN32 AND NOT UNIX)
  set(testPgSQLName "${testPgSQLName}.sh")
  set(testPgSQLHead "cd $$(dirname $$0)\n")
endif(WIN32 AND NOT UNIX)
set(testPgSQLBody "\"${PGREGRESS_BINARY}\" --inputdir=sql "
              "--bindir=\"${PGREGRESS_BINDIR}\" --user=postgres rdkit-91 "
              "props btree molgist bfpgist-91 sfpgist slfpgist fps reaction "
              "${avalonRegress} ${inchiRegress}\n")
file(WRITE ${testPgSQLName} ${testPgSQLHead} ${testPgSQLBody})
add_test(testPgSQL ${testPgSQLName})
#add_test(NAME testPgSQL
#              WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Code/PgSQL/rdkit
#              COMMAND ${PostgreSQL_ROOT}/bin/pg_regress --inputdir=sql
#              --bindir=${PostgreSQL_ROOT}/bin --user=postgres rdkit-91
#              props btree molgist bfpgist-91 sfpgist slfpgist fps reaction
#              ${avalonRegress} ${inchiRegress})
