add_custom_target(chemdraw_support ALL)
include(CMakePrintHelpers)

# The ChemDraw Library requires expat and expatpp to run.
# this has an include expat.h and expatpp.h that needs to be included
# Currently this is an OLD version of expatpp from source forge is used
# and included in this directory.
#
# For builds, we currently need a target_include_directories
# and will need to be fixed in the future
if(RDK_BUILD_CHEMDRAW_SUPPORT)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/chemdraw/CDXIO.h")
    set(RELEASE_NO "1.0.10")
    set(MD5 "ab7fc589de0c4574354434920f1b027c")
    downloadAndCheckMD5("https://codeload.github.com/Glysade/chemdraw/tar.gz/refs/tags/v${RELEASE_NO}"
      "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz" ${MD5})

    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
      ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-${RELEASE_NO}" "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw")
  endif()

  include(TestBigEndian)
  message("-- Looking for endianess")
  test_big_endian(WORDS_BIGENDIAN)

  # /* 1234 = LIL_ENDIAN, 4321 = BIGENDIAN */
  if(WORDS_BIGENDIAN)
    message("--   CHEDRAW BIGENDIAN PLATFORM")
    add_definitions("-DPLATFORM_BIGENDIAN")
  else(WORDS_BIGENDIAN)
    message("-    CHEMDRAW LITTLEENDIAN PLATFORM")
    add_definitions("-DPLATFORM_LITTLEENDIAN")
  endif(WORDS_BIGENDIAN)

  # we don't want to install expat, this is statically linked in to the ChemDraw lib
  # however, we don't want to install it so use the undocumented EXCLUDE_FROM_ALL
  add_subdirectory(chemdraw/expatpp EXCLUDE_FROM_ALL)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/expatpp/expatpp-code-r6-trunk/src_pp
    ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/expatpp/expatpp-code-r6-trunk/expat/lib)

  # it's way easier to use the RDKIT machinery to build and link so let's do that
  file(GLOB CHEMDRAW_SOURCE "chemdraw/chemdraw/*.cpp")
  rdkit_library(ChemDraw ${CHEMDRAW_SOURCE} SHARED)
  target_compile_definitions(ChemDraw PRIVATE CHEMDRAW_BUILD)
  target_link_libraries(ChemDraw PRIVATE expat)

  # export all the symbols for ChemDraw on MSVC
  if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR((NOT MSVC) AND WIN32))
    message("== ChemDraw exporting all symbols")
    set_target_properties(ChemDraw PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  endif()

  # On Windows, define MYLIB_EXPORTS when building the DLL
  if(WIN32)
    target_compile_definitions(ChemDraw
      PRIVATE CHEMDRAW_BUILD
    )
  endif()

  # On Linux/macOS, hide all symbols by default and expose only our API
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
    CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID STREQUAL "Emscripten")
    # Require -fvisibility=hidden
    set_target_properties(ChemDraw PROPERTIES
      CXX_VISIBILITY_PRESET hidden
      VISIBILITY_INLINES_HIDDEN ON
    )
  endif()

  # export all the symbols for ChemDraw on MSVC
  if((MSVC AND BUILD_SHARED_LIBS) OR((NOT MSVC) AND WIN32))
    set_target_properties(ChemDraw PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  endif()

  # Suppress warnings since we don't control the chedraw code and indicate
  # we have the EXPAT_CONFIG created
  if(MSVC)
    ADD_DEFINITIONS("-DTARGET_API_LIB -D_WINDOWS -DTARGET_OS_WIN32 -DHAVE_EXPAT_CONFIG_H")

    # we don't really control chemdraw source code, so suppress warnings
    target_compile_options(ChemDraw PRIVATE "/W0")
  else()
    ADD_DEFINITIONS("-DTARGET_API_LIB -D__linux -DHAVE_EXPAT_CONFIG_H")

    # we don't really control chemdraw source code, so suppress warnings
    target_compile_options(ChemDraw PRIVATE -w -Wno-unknown-pragmas -Wno-error)

    if(RDK_INSTALL_STATIC_LIBS)
      if(TARGET ChemDraw_static)
        target_compile_options(ChemDraw_static PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      endif()
    endif(RDK_INSTALL_STATIC_LIBS)
  endif()

  install(TARGETS ChemDraw DESTINATION ${RDKit_LibDir})

  # When building static libs, expat will not be included into ChemDraw,
  # so it needs to be installed so that it can be linked later on
  # when the libs are used.
  if(RDK_INSTALL_STATIC_LIBS OR NOT BUILD_SHARED_LIBS)
    install(TARGETS expat EXPORT ${RDKit_EXPORTED_TARGETS} COMPONENT dev)
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-comment -Wno-parentheses -Wno-logical-op-parentheses -Wno-pointer-bool-conversion -Wno-unused-value -Wno-unsequenced -Wno-constant-logical-operand")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wformat-overflow=0 -Wformat=0 -Wno-format-security")
  endif()

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
    SHARED LINK_LIBRARIES ChemDraw
    CIPLabeler ChemTransforms GraphMol RDGeneral Depictor SubstructMatch SmilesParse)

  rdkit_library(RDChemDrawReactionLib
    chemdrawreaction.cpp

    # ${EXPAT_SRC}
    SHARED LINK_LIBRARIES RDChemDrawLib ChemDraw
    CIPLabeler ChemTransforms ChemReactions GraphMol RDGeneral Depictor SubstructMatch SmilesParse)

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
    LINK_LIBRARIES RDChemDrawLib ChemDraw SubstructMatch ChemReactions
    FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawChiralCatchTest test-chiral.cpp
    LINK_LIBRARIES RDChemDrawLib ChemDraw SubstructMatch ChemReactions
    FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawReactionsCatchTest test-reactions.cpp
    LINK_LIBRARIES RDChemDrawReactionLib RDChemDrawLib ChemDraw SubstructMatch ChemReactions
    FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw3DCatchTest test_3d.cpp
    LINK_LIBRARIES RDChemDrawLib ChemDraw SubstructMatch ChemReactions
    FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw6KCatchTest test_6k.cpp
    LINK_LIBRARIES RDChemDrawLib ChemDraw SubstructMatch ChemReactions
    FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  if(RDK_BUILD_CPP_TESTS)
    if(MSVC)
      # The nanotubes blow up the smiles writer stack on MSVC so increase it
      set_target_properties(chemdraw6KCatchTest PROPERTIES LINK_FLAGS
        "/STACK:4194304")

    # this sets everything I think
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:4194304")
    else()
      # target_compile_options(chemdrawCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      # target_compile_options(chemdrawChiralCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      # target_compile_options(chemdrawReactionsCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      # target_compile_options(chemdraw3DCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      # target_compile_options(chemdraw6KCatchTest PRIVATE -w -Wno-unknown-pragmas -Wno-error)
    endif(MSVC)
  endif(RDK_BUILD_CPP_TESTS)

  if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)
endif(RDK_BUILD_CHEMDRAW_SUPPORT)
