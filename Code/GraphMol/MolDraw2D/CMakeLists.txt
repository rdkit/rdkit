
option(RDK_BUILD_QT_SUPPORT "build support for QT drawing" OFF )
option(RDK_BUILD_QT_DEMO "build the QT drawing demo" OFF )
option(RDK_BUILD_CAIRO_SUPPORT "build support for Cairo drawing" OFF )
option(RDK_BUILD_FREETYPE_SUPPORT "build support for FreeType font handling" ON)
option(RDK_INSTALL_COMIC_FONTS "download and install comic-neue to maximize the \"utility\" of the comic mode in MolDraw2D" ON)

if(RDK_BUILD_FREETYPE_SUPPORT AND RDK_INSTALL_COMIC_FONTS)
  set(needDownload "TRUE")
  if(EXISTS "${RDKit_DataDir}/Fonts/ComicNeue-Regular.ttf")
    message("${RDKit_DataDir}/Fonts/ComicNeue-Regular.ttf already there, skipping the download")
    set(needDownload "FALSE")
  endif()
  if(needDownload)
    set(MD5Sum "23ed3f833c1ae0adb141a26b4a30d73e")
    downloadAndCheckMD5("https://fonts.google.com/download?family=Comic%20Neue"
          "${CMAKE_CURRENT_SOURCE_DIR}/Comic_Neue.zip"
          ${MD5Sum})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar x
      ${CMAKE_CURRENT_SOURCE_DIR}/Comic_Neue.zip --format=zip
      WORKING_DIRECTORY ${RDKit_DataDir}/Fonts)
  endif(needDownload)
endif()

rdkit_headers(MolDraw2D.h
  MolDraw2DSVG.h
  MolDraw2Dwx.h
  MolDraw2DJS.h
  MolDraw2DUtils.h
  DrawText.h
  DrawTextSVG.h
  DrawTextJS.h
  DEST GraphMol/MolDraw2D
)

rdkit_library(MolDraw2D MolDraw2D.cpp MolDraw2DSVG.cpp
  MolDraw2DDetails.cpp MolDraw2DUtils.cpp DrawText.cpp
  DrawTextSVG.cpp 
  MolDraw2DJS.cpp DrawTextJS.cpp
  LINK_LIBRARIES
  ChemReactions Depictor MolTransforms 
  SubstructMatch Subgraphs GraphMol EigenSolvers )
target_compile_definitions(MolDraw2D PRIVATE RDKIT_MOLDRAW2D_BUILD)

# if(RDK_BUILD_QT_DEMO && !RDK_BUILD_QT_SUPPORT) error
if(RDK_BUILD_QT_SUPPORT)
  add_subdirectory(Qt)
endif(RDK_BUILD_QT_SUPPORT)

if(RDK_BUILD_CAIRO_SUPPORT)
  find_package(Cairo REQUIRED)
  target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
  target_link_libraries(MolDraw2D PUBLIC Cairo::Cairo)
  target_sources(MolDraw2D PRIVATE MolDraw2DCairo.cpp DrawTextCairo.cpp)
  if(RDK_INSTALL_STATIC_LIBS)
     target_compile_definitions(MolDraw2D_static PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
     target_link_libraries(MolDraw2D_static PUBLIC Cairo::Cairo)
     target_sources(MolDraw2D_static PRIVATE MolDraw2DCairo.cpp DrawTextCairo.cpp)
  endif()
  rdkit_headers(MolDraw2DCairo.h DrawTextCairo.h DEST GraphMol/MolDraw2D)
endif(RDK_BUILD_CAIRO_SUPPORT)

if(RDK_BUILD_FREETYPE_SUPPORT)
  target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
  target_sources(MolDraw2D PRIVATE DrawText.cpp
          DrawTextFT.cpp DrawTextFTSVG.cpp DrawTextFTJS.cpp)
  if(RDK_INSTALL_STATIC_LIBS)
    target_compile_definitions(MolDraw2D_static PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
    target_sources(MolDraw2D_static PRIVATE DrawText.cpp
          DrawTextFT.cpp DrawTextFTSVG.cpp DrawTextFTJS.cpp)
  endif()
  if( ${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
      set(USE_FLAGS "-s USE_FREETYPE=1")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${USE_FLAGS}")
      target_link_libraries(MolDraw2D PUBLIC ${FREETYPE_LIBRARIES})
      if(RDK_INSTALL_STATIC_LIBS)
        target_link_libraries(MolDraw2D_static PUBLIC ${FREETYPE_LIBRARIES})
      endif()
  else()
      find_package(Freetype REQUIRED)
      target_link_libraries(MolDraw2D PUBLIC Freetype::Freetype)
      if(RDK_INSTALL_STATIC_LIBS)
        target_link_libraries(MolDraw2D_static PUBLIC Freetype::Freetype)
      endif()
  endif()
  rdkit_headers(DrawText.h DrawTextFT.h
          DrawTextFTSVG.h DrawTextFTJS.h DEST GraphMol/MolDraw2D)
  if(RDK_BUILD_CAIRO_SUPPORT)
    target_sources(MolDraw2D PRIVATE DrawTextFTCairo.cpp)
    if(RDK_INSTALL_STATIC_LIBS)
      target_sources(MolDraw2D_static PRIVATE DrawTextFTCairo.cpp)
    endif()
    rdkit_headers(DrawTextFTCairo.h DEST GraphMol/MolDraw2D)
  endif()
endif(RDK_BUILD_FREETYPE_SUPPORT)

rdkit_test(moldraw2DTest1 test1.cpp LINK_LIBRARIES
  MolDraw2D )

rdkit_catch_test(moldraw2DTestCatch ../catch_main.cpp catch_tests.cpp LINK_LIBRARIES
  MolDraw2D CIPLabeler )

rdkit_test(moldraw2DRxnTest1 rxn_test1.cpp LINK_LIBRARIES
  MolDraw2D )

  if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
