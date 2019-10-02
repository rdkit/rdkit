add_custom_target(inchi_support ALL)

if(RDK_BUILD_INCHI_SUPPORT)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/INCHI_BASE/src/ichican2.c)
    add_definitions(-DRDKIT_INCHI_BUILD)
    rdkit_library(Inchi
      src/INCHI_BASE/src/ichi_bns.c
      src/INCHI_BASE/src/ichi_io.c
      src/INCHI_BASE/src/ichican2.c
      src/INCHI_BASE/src/ichicano.c
      src/INCHI_BASE/src/ichicans.c
      src/INCHI_BASE/src/ichierr.c
      src/INCHI_BASE/src/ichiisot.c
      src/INCHI_BASE/src/ichimak2.c
      src/INCHI_BASE/src/ichimake.c
      src/INCHI_BASE/src/ichimap1.c
      src/INCHI_BASE/src/ichimap2.c
      src/INCHI_BASE/src/ichimap4.c
      src/INCHI_BASE/src/ichinorm.c
      src/INCHI_BASE/src/ichiparm.c
      src/INCHI_BASE/src/ichiprt1.c
      src/INCHI_BASE/src/ichiprt2.c
      src/INCHI_BASE/src/ichiprt3.c
      src/INCHI_BASE/src/ichiqueu.c
      src/INCHI_BASE/src/ichiread.c
      src/INCHI_BASE/src/ichiring.c
      src/INCHI_BASE/src/ichirvr1.c
      src/INCHI_BASE/src/ichirvr2.c
      src/INCHI_BASE/src/ichirvr3.c
      src/INCHI_BASE/src/ichirvr4.c
      src/INCHI_BASE/src/ichirvr5.c
      src/INCHI_BASE/src/ichirvr6.c
      src/INCHI_BASE/src/ichirvr7.c
      src/INCHI_BASE/src/ichisort.c
      src/INCHI_BASE/src/ichister.c
      src/INCHI_BASE/src/ichitaut.c
      src/INCHI_BASE/src/ikey_base26.c
      src/INCHI_BASE/src/ikey_dll.c
      src/INCHI_BASE/src/mol2atom.c
      src/INCHI_BASE/src/mol_fmt1.c
      src/INCHI_BASE/src/mol_fmt2.c
      src/INCHI_BASE/src/mol_fmt3.c
      src/INCHI_BASE/src/mol_fmt4.c
      src/INCHI_BASE/src/readinch.c
      src/INCHI_BASE/src/runichi.c
      src/INCHI_BASE/src/runichi2.c
      src/INCHI_BASE/src/runichi3.c
      src/INCHI_BASE/src/runichi4.c
      src/INCHI_BASE/src/sha2.c
      src/INCHI_BASE/src/strutil.c
      src/INCHI_BASE/src/util.c
      src/INCHI_API/libinchi/src/ichilnct.c
      src/INCHI_API/libinchi/src/inchi_dll.c
      src/INCHI_API/libinchi/src/inchi_dll_a.c
      src/INCHI_API/libinchi/src/inchi_dll_a2.c
      src/INCHI_API/libinchi/src/inchi_dll_b.c
      src/INCHI_API/libinchi/src/inchi_dll_main.c
      src/INCHI_API/libinchi/src/ixa/ixa_builder.c
      src/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c
      src/INCHI_API/libinchi/src/ixa/ixa_mol.c
      src/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c
      src/INCHI_API/libinchi/src/ixa/ixa_read_mol.c
      src/INCHI_API/libinchi/src/ixa/ixa_status.c SHARED)
      if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
        set_target_properties(Inchi PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
      endif()
    install(TARGETS Inchi DESTINATION ${RDKit_LibDir})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/INCHI_BASE/src/)
    ADD_DEFINITIONS("-DTARGET_API_LIB")
  endif()

     if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-parentheses -Wno-logical-op-parentheses")
     endif()
  add_definitions(-DRDKIT_RDINCHILIB_BUILD)
  rdkit_library(RDInchiLib inchi.cpp SHARED LINK_LIBRARIES ${INCHI_LIBRARIES}
                GraphMol RDGeneral Depictor SubstructMatch SmilesParse
                ${RDKit_THREAD_LIBS})
  rdkit_headers(inchi.h DEST GraphMol)
  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  add_definitions(-DRDK_BUILD_INCHI_SUPPORT)

  rdkit_test(testInchi test.cpp
    LINK_LIBRARIES
    ${INCHI_LIBRARIES}
    RDInchiLib ${INCHI_LIBRARIES}
    FileParsers SmilesParse Descriptors Depictor SubstructMatch GraphMol
    RDGeneral DataStructs RDGeneral RDGeometryLib
    ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_INCHI_SUPPORT)

if(RDK_BUILD_PYTHON_WRAPPERS)
  # we need an inchi.py... copy in the appropriate one
  if(RDK_BUILD_INCHI_SUPPORT)
   add_custom_command(TARGET inchi_support
                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/inchi.py
                     ${RDK_PYTHON_BINARY_DIR}/rdkit/Chem/inchi.py)
  else(RDK_BUILD_INCHI_SUPPORT)
    add_custom_command(TARGET inchi_support
                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/noinchi.py
                     ${RDK_PYTHON_BINARY_DIR}/rdkit/Chem/inchi.py)
  endif(RDK_BUILD_INCHI_SUPPORT)
  get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
          "${extra_clean_files};${RDK_PYTHON_BINARY_DIR}/rdkit/Chem/inchi.py")
endif(RDK_BUILD_PYTHON_WRAPPERS)
