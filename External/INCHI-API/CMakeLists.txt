add_custom_target(inchi_support ALL)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)
	add_custom_command(
		TARGET inchi_support
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/inchi.py ${CMAKE_BINARY_DIR}/rdkit/Chem/inchi.py)
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

	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
	rdkit_library(RDInchiLib inchi.cpp SHARED LINK_LIBRARIES Inchi GraphMol
		      RDGeneral Depictor SubstructMatch)

if(RDK_BUILD_PYTHON_WRAPPERS)
	add_subdirectory(Wrap)
endif(RDK_BUILD_PYTHON_WRAPPERS)

        set(RDK_BUILDING_INCHI_SUPPORT ON)
        add_definitions(-DBUILD_INCHI_SUPPORT)
else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)

	message(WARNING  "** NO INCHI SOFTWARE FOUND\n"
	       	"InChI support will be turned off. If you want to add InChI support, please follow instructions inside External/INCHI-API/README to download a copy of InChI software and then rerun cmake.")

	add_custom_command(
		TARGET inchi_support
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/noinchi.py ${CMAKE_BINARY_DIR}/rdkit/Chem/inchi.py)
        set(RDK_BUILDING_INCHI_SUPPORT OFF)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/ichican2.c)

if(RDK_BUILD_PYTHON_WRAPPERS)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inchi.py DESTINATION ${RDKit_PythonDir}/Chem)
endif(RDK_BUILD_PYTHON_WRAPPERS)
