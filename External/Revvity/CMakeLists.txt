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
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/CDXIO.h)
    # Build expat
    # set(BUILD_shared ON)
    # add_subdirectory(expatpp/expatpp-code-r6-trunk/expat)
    # Fix bug in the expat 2.1.0 CMakeFile on windows
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
    else()
      ADD_DEFINITIONS("-DTARGET_API_LIB -D__linux -DHAVE_EXPAT_CONFIG_H")
    endif()
      
    install(TARGETS ChemDraw DESTINATION ${RDKit_LibDir})

    # get the big endiand flags
    test_big_endian(WORDS_BIGENDIAN)
    #/* 1234 = LIL_ENDIAN, 4321 = BIGENDIAN */
    if(WORDS_BIGENDIAN)
          add_compile_definitions("PLATFORM_BIGENDIAN")	
    else(WORDS_BIGENDIAN)
      add_compile_definitions("PLATFORM_LITTLEENDIAN")
    endif(WORDS_BIGENDIAN)

  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-comment -Wno-parentheses -Wno-logical-op-parentheses -Wno-pointer-bool-conversion -Wno-unused-value -Wno-unsequenced -Wno-constant-logical-operand")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wformat-overflow=0 -Wformat=0 -Wno-format-security")
  endif()

  # finally the RDKit library!
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
          CIPLabeler ChemTransforms ChemReactions GraphMol RDGeneral Depictor SubstructMatch SmilesParse )

  if(MSVC)
      target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD XML_USE_MSC_EXTENSIONS)
  else()
    target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD)
  endif()

  rdkit_headers(chemdraw.h DEST GraphMol)
  if(RDK_BUILD_PYTHON_WRAPPERS)
    #add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  add_definitions(-DRDK_BUILD_CHEMDRAW_SUPPORT)

  # all the tests
  rdkit_catch_test(chemdrawCatchTest test.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawChiralCatchTest test-chiral.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdrawReactionsCatchTest test-reactions.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw3DCatchTest test_3d.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw6KCatchTest test_6k.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch ChemReactions FileParsers SmilesParse CIPLabeler ChemTransforms GraphMol)

  if(MSVC)
    # The nanotubes blow up the smiles writer stack on MSVC so increase it
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:4194304")
  endif(MSVC)

endif(RDK_BUILD_CHEMDRAW_SUPPORT)

