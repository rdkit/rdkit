add_custom_target(chemdraw_support ALL)
include(CMakePrintHelpers)

if(RDK_BUILD_CHEMDRAW_SUPPORT)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/CDXIO.h)
    add_subdirectory(expatpp/expatpp-code-r6-trunk/expat)
    rdkit_library(ChemDraw
        src/Base64Converter.cpp
        src/CDXAnnotation.cpp
        src/CDXArrow.cpp
        src/CDXAttachedDataHelper.cpp
        src/CDXBandMarker.cpp
        src/CDXBioDraw.cpp
        src/CDXBond.cpp
        src/CDXBorder.cpp
        src/CDXBracketedGroup.cpp
        src/CDXChemicalProperty.cpp
        src/CDXColorTable.cpp
        src/CDXColoredMolecularArea.cpp
        src/CDXConstraint.cpp
        src/CDXCurve.cpp
        src/CDXDocument.cpp
        src/CDXDocumentProperty.cpp
        src/CDXDocumentPropertyCollection.cpp
        src/CDXEmbeddedObject.cpp
        src/CDXFontTable.cpp
        src/CDXFragment.cpp
        src/CDXGEPBand.cpp
        src/CDXGEPPlate.cpp
        src/CDXGeometry.cpp
        src/CDXGraphic.cpp
        src/CDXGraphicObject.cpp
        src/CDXGroup.cpp
        src/CDXGroupOrFragment.cpp
        src/CDXIO.cpp
        src/CDXIndicatorText.cpp
        src/CDXMLParser.cpp
        src/CDXNamedAlternativeGroup.cpp
        src/CDXNode.cpp
        src/CDXObject.cpp
        src/CDXObjectTag.cpp
        src/CDXPage.cpp
        src/CDXPlasmidMap.cpp
        src/CDXPlateBase.cpp
        src/CDXPlateItemBase.cpp
        src/CDXRLogic.cpp
        src/CDXReactionStep.cpp
        src/CDXSpectrum.cpp
        src/CDXSplitter.cpp
        src/CDXStdObjects.cpp
        src/CDXStoichiometryGrid.cpp
        src/CDXTLCLane.cpp
        src/CDXTLCPlate.cpp
        src/CDXTLCSpot.cpp
        src/CDXTable.cpp
        src/CDXTemplateGrid.cpp
        src/CDXTest.cpp
        src/CDXText.cpp
        src/CDXUnicode.cpp
        src/CDXUtils.cpp
        src/CDXUtils_Win.cpp
        src/XMLDoc.cpp
        src/XMLParser.cpp
	expatpp/expatpp-code-r6-trunk/expat/lib/xmlparse.c
	expatpp/expatpp-code-r6-trunk/expat/lib/xmlrole.c
	expatpp/expatpp-code-r6-trunk/expat/lib/xmltok.c
	expatpp/expatpp-code-r6-trunk/expat/lib/xmltok_impl.c
	expatpp/expatpp-code-r6-trunk/expat/lib/xmltok_ns.c
	expatpp/expatpp-code-r6-trunk/src_pp/expatpp.cpp
	)
	
    target_compile_definitions(ChemDraw PRIVATE RDKIT_CHEMDRAW_BUILD)

    if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
        set_target_properties(ChemDraw PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    endif()
      
    install(TARGETS ChemDraw DESTINATION ${RDKit_LibDir})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
                        ${CMAKE_CURRENT_SOURCE_DIR}/expatpp/expatpp-code-r6-trunk/src_pp)

    test_big_endian(WORDS_BIGENDIAN)
    #/* 1234 = LIL_ENDIAN, 4321 = BIGENDIAN */
    if(WORDS_BIGENDIAN)
          add_compile_definitions("PLATFORM_BIGENDIAN")	
    else(WORDS_BIGENDIAN)
      add_compile_definitions("PLATFORM_LITTLEENDIAN")
    endif(WORDS_BIGENDIAN)

    ADD_DEFINITIONS("-DTARGET_API_LIB -D__linux" "-DHAVE_EXPAT_CONFIG_H")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-comment -Wno-parentheses -Wno-logical-op-parentheses -Wno-pointer-bool-conversion -Wno-unused-value -Wno-unsequenced -Wno-constant-logical-operand")
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wformat-overflow=0 -Wformat=0 -Wno-format-security")
  endif()

  rdkit_library(RDChemDrawLib
      bond.cpp
      chemdraw.cpp
      fragment.cpp
      node.cpp
      reaction.cpp
      utils.cpp
      SHARED LINK_LIBRARIES ChemDraw
          ChemTransforms GraphMol RDGeneral Depictor SubstructMatch SmilesParse )
  target_compile_definitions(RDChemDrawLib PRIVATE RDKIT_RDCHEMDRAWLIB_BUILD)

  rdkit_headers(chemdraw.h DEST GraphMol)
  if(RDK_BUILD_PYTHON_WRAPPERS)
    #add_subdirectory(Wrap)
  endif(RDK_BUILD_PYTHON_WRAPPERS)

  add_definitions(-DRDK_BUILD_CHEMDRAW_SUPPORT)

  rdkit_catch_test(chemdrawCatchTest test.cpp
        LINK_LIBRARIES RDChemDrawLib SubstructMatch FileParsers SmilesParse ChemTransforms GraphMol)

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
