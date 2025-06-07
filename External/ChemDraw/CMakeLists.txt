add_custom_target(chemdraw_support ALL)
include(CMakePrintHelpers)

# The ChemDraw Library requires expat and expatpp to run.
#  this has an include expat.h and expatpp.h that needs to be included
#  Currently this is an OLD version of expatpp from source forge is used
#  and included in this directory.
# 
#  For builds, we currently need a target_include_directories
#  and will need to be fixed in the future

if(RDK_BUILD_CHEMDRAW_SUPPORT)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/chemdraw/CDXIO.h" )
    set(RELEASE_NO "1.0.7")
    set(MD5 "a7d5c49a21384278698946fb1051b474")
    downloadAndCheckMD5("https://codeload.github.com/Glysade/chemdraw/tar.gz/refs/tags/v${RELEASE_NO}"
              "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz" ${MD5})
	      
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-${RELEASE_NO}" "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw")
   endif()
   
   add_subdirectory(chemdraw)
   include_directories(chemdraw)
   
     set(RDChemDrawLib_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        CACHE STRING "RDChemDrawLib Include File" FORCE)

  rdkit_library(RDChemDrawLib
      bond.cpp
      bracket.cpp
      chemdraw.cpp
      fragment.cpp
      node.cpp
      reaction.cpp
      utils.cpp
      writer.cpp
      # ${EXPAT_SRC}
      SHARED LINK_LIBRARIES ChemDraw expat
          CIPLabeler ChemTransforms GraphMol RDGeneral Depictor SubstructMatch SmilesParse )

  rdkit_library(RDChemDrawReactionLib
      chemdrawreaction.cpp
      # ${EXPAT_SRC}
      SHARED LINK_LIBRARIES RDChemDrawLib ChemDraw expat
          CIPLabeler ChemTransforms ChemReactions GraphMol RDGeneral Depictor SubstructMatch SmilesParse )

  if(MSVC)
      target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD XML_USE_MSC_EXTENSIONS)
      target_compile_definitions(RDChemDrawReactionLib PRIVATE RDKIT_RDCHEMDRAWREACTIONLIB_BUILD
                                                               XML_USE_MSC_EXTENSIONS)
  else()
    target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD)
    target_compile_definitions(RDChemDrawReactionLib PRIVATE RDKIT_RDCHEMDRAWREACTIONLIB_BUILD)
  endif()

  install(TARGETS RDChemDrawLib DESTINATION ${RDKit_LibDir})
  install(TARGETS RDChemDrawReactionLib DESTINATION ${RDKit_LibDir})
  set(RDK_CHEMDRAW_LIBS RDChemDrawLib CACHE STRING "the external libraries" FORCE)
  set(RDK_CHEMDRAWREACTION_LIBS RDChemDrawReactinoLib CACHE STRING "the external libraries" FORCE)  

  rdkit_headers(chemdraw.h DEST GraphMol)
  rdkit_headers(chemdrawreaction.h DEST GraphMol)

  # all the tests
  rdkit_catch_test(chemdrawCatchTest test.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions
	FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawChiralCatchTest test-chiral.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions
	FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawReactionsCatchTest test-reactions.cpp
        LINK_LIBRARIES RDChemDrawReactionLib RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions
	FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw3DCatchTest test_3d.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions
	FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw6KCatchTest test_6k.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions
	FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  if(RDK_BUILD_CPP_TESTS)
    if(MSVC)
    # The nanotubes blow up the smiles writer stack on MSVC so increase it
    set_target_properties(chemdraw6KCatchTest PROPERTIES LINK_FLAGS
                   "/STACK:4194304")

    # this sets everything I think
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:4194304")
    else()
      #target_compile_options(chemdrawCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      #target_compile_options(chemdrawChiralCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      #target_compile_options(chemdrawReactionsCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      #target_compile_options(chemdraw3DCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      #target_compile_options(chemdraw6KCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)  
    endif(MSVC)
  endif(RDK_BUILD_CPP_TESTS)

  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)


endif(RDK_BUILD_CHEMDRAW_SUPPORT)

