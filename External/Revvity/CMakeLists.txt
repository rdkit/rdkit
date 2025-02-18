add_custom_target(chemdraw_support ALL)
include(CMakePrintHelpers)

if(RDK_BUILD_CHEMDRAW_SUPPORT)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/CDXIO.h)
    # Build expat
    # set(BUILD_shared ON)
    # add_subdirectory(expatpp/expatpp-code-r6-trunk/expat)
    # Fix bug in the expat 2.1.0 CMakeFile on windows
    add_subdirectory(expatpp)
    
    target_include_directories(expat PRIVATE
                        ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/
                        ${CMAKE_CURRENT_SOURCE_DIR}/expat/
			${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                        ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp
			)


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
    target_include_directories(ChemDraw PRIVATE
                        ${CMAKE_CURRENT_SOURCE_DIR}/chemdraw/
                        ${CMAKE_CURRENT_SOURCE_DIR}/expat/
			${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                        ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp
			)

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
          ChemTransforms GraphMol RDGeneral Depictor SubstructMatch SmilesParse )

  if(MSVC)
      target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD XML_USE_MSC_EXTENSIONS)
  else()
    target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD)
  endif()

  target_include_directories(RDChemDrawLib PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/expat/
		     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp/)


  rdkit_headers(chemdraw.h DEST GraphMol)
  if(RDK_BUILD_PYTHON_WRAPPERS)
    #add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  add_definitions(-DRDK_BUILD_CHEMDRAW_SUPPORT)

  rdkit_catch_test(chemdrawCatchTest test.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch FileParsers SmilesParse ChemTransforms GraphMol)
      
  rdkit_catch_test(chemdraw3DCatchTest test_3d.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch FileParsers SmilesParse ChemTransforms GraphMol)

  rdkit_catch_test(chemdraw6KCatchTest test_6k.cpp
        LINK_LIBRARIES RDChemDrawLib ChemDraw expat SubstructMatch FileParsers SmilesParse ChemTransforms GraphMol)

  target_include_directories(chemdrawCatchTest PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/expat/
		     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp/)

  target_include_directories(chemdraw3DCatchTest PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/expat/
		     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp/)

  target_include_directories(chemdraw6KCatchTest PRIVATE
                     ${CMAKE_CURRENT_SOURCE_DIR}/expat/
		     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/expat/lib
                     ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp/)
endif(RDK_BUILD_CHEMDRAW_SUPPORT)

if(RDK_BUILD_PYTHON_WRAPPERS)
  # we need an inchi.py... copy in the appropriate one
#  if(RDK_BUILD_CHEMDRAW_SUPPORT)
#   add_custom_command(TARGET inchi_support
#                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/inchi.py
#                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
#  else(RDK_BUILD_CHEMDRAW_SUPPORT)
#    add_custom_command(TARGET inchi_support
#                  COMMAND ${C*MAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/noinchi.py
#                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
#  endif(RDK_BUILD_CHEMDRAW_SUPPORT)
#  get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
#  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
#          "${extra_clean_files};${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py")
endif(RDK_BUILD_PYTHON_WRAPPERS)
