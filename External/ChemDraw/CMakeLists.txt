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
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/CDXIO.h" )
    set(RELEASE_NO "1.0.3")
    set(MD5 "b3a79d047b4aab69cd9935be1f038e1c")
    downloadAndCheckMD5("https://codeload.github.com/Glysade/chemdraw/tar.gz/refs/tags/v${RELEASE_NO}"
              "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz" ${MD5})
	      
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-v${RELEASE_NO}.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-${RELEASE_NO}/chemdraw" "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw")
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/chemdraw-${RELEASE_NO}/expatpp" "${CMAKE_CURRENT_SOURCE_DIR}/expatpp")
   endif()
   
    # Build expat
    add_subdirectory(expatpp)

    #
    # ChemDraw library for reading/writing documents
    rdkit_library(ChemDraw
        chemdraw/Base64Converter.cpp
        chemdraw/CDXAnnotation.cpp
        chemdraw/CDXArrow.cpp
        chemdraw/CDXAttachedDataHelper.cpp
        chemdraw/CDXBandMarker.cpp
        chemdraw/CDXBioDraw.cpp
        chemdraw/CDXBond.cpp
        chemdraw/CDXBorder.cpp
        chemdraw/CDXBracketedGroup.cpp
        chemdraw/CDXChemicalProperty.cpp
        chemdraw/CDXColorTable.cpp
        chemdraw/CDXColoredMolecularArea.cpp
        chemdraw/CDXConstraint.cpp
        chemdraw/CDXCurve.cpp
        chemdraw/CDXDocument.cpp
        chemdraw/CDXDocumentProperty.cpp
        chemdraw/CDXDocumentPropertyCollection.cpp
        chemdraw/CDXEmbeddedObject.cpp
        chemdraw/CDXFontTable.cpp
        chemdraw/CDXFragment.cpp
        chemdraw/CDXGEPBand.cpp
        chemdraw/CDXGEPPlate.cpp
        chemdraw/CDXGeometry.cpp
        chemdraw/CDXGraphic.cpp
        chemdraw/CDXGraphicObject.cpp
        chemdraw/CDXGroup.cpp
        chemdraw/CDXGroupOrFragment.cpp
        chemdraw/CDXIO.cpp
        chemdraw/CDXIndicatorText.cpp
        chemdraw/CDXMLParser.cpp
        chemdraw/CDXNamedAlternativeGroup.cpp
        chemdraw/CDXNode.cpp
        chemdraw/CDXObject.cpp
        chemdraw/CDXObjectTag.cpp
        chemdraw/CDXPage.cpp
        chemdraw/CDXPlasmidMap.cpp
        chemdraw/CDXPlateBase.cpp
        chemdraw/CDXPlateItemBase.cpp
        chemdraw/CDXRLogic.cpp
        chemdraw/CDXReactionStep.cpp
        chemdraw/CDXSpectrum.cpp
        chemdraw/CDXSplitter.cpp
        chemdraw/CDXStdObjects.cpp
        chemdraw/CDXStoichiometryGrid.cpp
        chemdraw/CDXTLCLane.cpp
        chemdraw/CDXTLCPlate.cpp
        chemdraw/CDXTLCSpot.cpp
        chemdraw/CDXTable.cpp
        chemdraw/CDXTemplateGrid.cpp
        chemdraw/CDXTest.cpp
        chemdraw/CDXText.cpp
        chemdraw/CDXUnicode.cpp
        chemdraw/CDXUtils.cpp
        chemdraw/CDXUtils_Win.cpp
        chemdraw/XMLDoc.cpp
        chemdraw/XMLParser.cpp
	# ${EXPAT_SRC}
	SHARED LINK_LIBRARIES expat
	)

    target_compile_definitions(ChemDraw PRIVATE RDKIT_CHEMDRAW_BUILD)
    

    # export all the symbols for ChemDraw on MSVC
    if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
        message("== ChemDraw exporting all symbols")
        set_target_properties(ChemDraw PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    endif()
    
    if(MSVC)
      ADD_DEFINITIONS("-DTARGET_API_LIB -D_WINDOWS -DTARGET_OS_WIN32 -DHAVE_EXPAT_CONFIG_H")
      # we don't really control chemdraw source code, so suppress warnings	
      target_compile_options(ChemDraw PRIVATE "/W0")
    else()
      ADD_DEFINITIONS("-DTARGET_API_LIB -D__linux -DHAVE_EXPAT_CONFIG_H")
      # we don't really control chemdraw source code, so suppress warnings
      target_compile_options(ChemDraw PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      if(RDK_INSTALL_STATIC_LIBS)
            target_compile_options(ChemDraw_static PRIVATE -w -Wno-unknown-pragmas -Wno-error)
      endif(RDK_INSTALL_STATIC_LIBS)
    endif()
      
    install(TARGETS ChemDraw DESTINATION ${RDKit_LibDir})



  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-comment -Wno-parentheses -Wno-logical-op-parentheses -Wno-pointer-bool-conversion -Wno-unused-value -Wno-unsequenced -Wno-constant-logical-operand")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wformat-overflow=0 -Wformat=0 -Wno-format-security")
  endif()

  #######################################################################################
  # finally the RDKit library!
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

