add_custom_target(inchi_support ALL)

if(RDK_BUILD_INCHI_SUPPORT)
  # search system InChI
  find_package(Inchi)

  # check whether we have custom InChI source
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)
    message(STATUS  "Found InChI software locally")
    if(INCHI_FOUND)
      message(WARNING "** Local InChI software takes precedence when both system InChI and local InChI are found")
    endif(INCHI_FOUND)
    rdkit_library(Inchi src/ichican2.c       src/ichicano.c 
                       src/ichicans.c       src/ichiisot.c 
                       src/ichilnct.c       src/ichimak2.c
                       src/ichimake.c       src/ichimap1.c 
                       src/ichimap2.c       src/ichimap4.c 
                       src/ichinorm.c       src/ichiparm.c 
                       src/ichiprt1.c       src/ichiprt2.c 
                       src/ichiprt3.c       src/ichiqueu.c 
                       src/ichiring.c       src/ichisort.c 
                       src/ichister.c       src/ichitaut.c 
                       src/ichi_bns.c       src/inchi_dll.c
                       src/ichiread.c       src/ichirvr1.c 
                       src/ichirvr2.c       src/ichirvr3.c 
                       src/ichirvr4.c       src/ichirvr5.c
                       src/ichirvr6.c       src/ichirvr7.c 
                       src/inchi_dll_main.c src/runichi.c  
                       src/inchi_dll_a.c    src/inchi_dll_a2.c  
                       src/ikey_dll.c       src/ikey_base26.c  
                       src/ichi_io.c  
                       src/sha2.c           src/strutil.c  
                       src/util.c SHARED)
    install(TARGETS Inchi DESTINATION ${RDKit_LibDir})
    set(INCHI_LIBRARIES Inchi)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

  else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)
    if (INCHI_FOUND)
      include_directories(${INCHI_INCLUDE_DIR})
    else (INCHI_FOUND)
      # system InChI is missing
      message(WARNING  "** NO INCHI SOFTWARE FOUND\n"
              "InChI support will be turned off. If you want to add InChI support, please follow the instructions inside $RDBASE/External/INCHI-API/README to download a copy of InChI software and then rerun cmake.")
      set(RDK_BUILD_INCHI_SUPPORT OFF)
    endif(INCHI_FOUND)
  endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)

endif(RDK_BUILD_INCHI_SUPPORT)

# check whether InChI support is ON after the search
if(RDK_BUILD_INCHI_SUPPORT)
  rdkit_library(RDInchiLib inchi.cpp SHARED LINK_LIBRARIES ${INCHI_LIBRARIES} 
                GraphMol RDGeneral Depictor SubstructMatch SmilesParse 
                ${Boost_THREAD_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})           
  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  add_definitions(-DBUILD_INCHI_SUPPORT)

rdkit_test(testInchi test.cpp 
LINK_LIBRARIES 
${INCHI_LIBRARIES} 
   RDInchiLib ${INCHI_LIBRARIES} 
   FileParsers SmilesParse Descriptors Depictor SubstructMatch GraphMol 
   RDGeneral DataStructs RDGeneral RDGeometryLib 
   ${Boost_THREAD_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
endif(RDK_BUILD_INCHI_SUPPORT)

if(RDK_BUILD_PYTHON_WRAPPERS)
  # we need an inchi.py... copy in the appropriate one
  if(RDK_BUILD_INCHI_SUPPORT)
   add_custom_command(TARGET inchi_support
                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/inchi.py
                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
  else(RDK_BUILD_INCHI_SUPPORT)
    add_custom_command(TARGET inchi_support
                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/noinchi.py
                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
  endif(RDK_BUILD_INCHI_SUPPORT)
  get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
          "${extra_clean_files};${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py")
endif(RDK_BUILD_PYTHON_WRAPPERS)

