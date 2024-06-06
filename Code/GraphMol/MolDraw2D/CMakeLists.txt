
option(RDK_BUILD_QT_SUPPORT "build support for QT drawing" OFF)
option(RDK_USE_QT6 "Use Qt6 instead of Qt5" OFF)
option(RDK_BUILD_QT_DEMO "build the QT drawing demo" OFF)
option(RDK_BUILD_CAIRO_SUPPORT "build support for Cairo drawing" OFF)
option(RDK_BUILD_FREETYPE_SUPPORT "build support for FreeType font handling" ON)
option(RDK_INSTALL_COMIC_FONTS "download and install comic-neue to maximize the \"utility\" of the comic mode in MolDraw2D" ON)

if (RDK_BUILD_MINIMAL_LIB AND RDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS)
    add_definitions(-DRDK_MINIMAL_LIB_SUPPORT_LEGACY_BROWSERS)
endif ()

if (RDK_BUILD_FREETYPE_SUPPORT AND RDK_INSTALL_COMIC_FONTS)
    set(needDownload "TRUE")
    if (EXISTS "${RDKit_DataDir}/Fonts/ComicNeue-Regular.ttf")
        message("${RDKit_DataDir}/Fonts/ComicNeue-Regular.ttf already there, skipping the download")
        set(needDownload "FALSE")
    endif ()
    if (needDownload)
        # we started having problems with constantly changing MD5s on the zip file,
        # so now we just check the MD5 of the target file
        downloadAndCheckMD5("https://github.com/google/fonts/raw/main/ofl/comicneue/ComicNeue-Regular.ttf"
                "${RDKit_DataDir}/Fonts/ComicNeue-Regular.ttf"
                "fc1eac54b325542d4c133732658f823b")
        downloadAndCheckMD5("https://github.com/google/fonts/raw/main/ofl/comicneue/OFL.txt"
                "${RDKit_DataDir}/Fonts/OFL.txt"
                "")
    endif (needDownload)

endif ()

rdkit_headers(MolDraw2D.h
        MolDraw2DSVG.h
        MolDraw2Dwx.h
        MolDraw2DJS.h
        MolDraw2DUtils.h
        MolDraw2DHelpers.h
        DEST GraphMol/MolDraw2D
)

rdkit_library(MolDraw2D MolDraw2D.cpp MolDraw2DSVG.cpp
        MolDraw2DDetails.cpp MolDraw2DUtils.cpp AtomSymbol.cpp
        DrawAnnotation.cpp DrawMol.cpp DrawMolMCH.cpp DrawMolMCHCircleAndLine.cpp
        DrawMolMCHLasso.cpp
        DrawShape.cpp DrawText.cpp DrawTextNotFT.cpp DrawTextSVG.cpp
        MolDraw2DJS.cpp DrawTextJS.cpp
        LINK_LIBRARIES
        ChemReactions Depictor MolTransforms
        SubstructMatch Subgraphs GraphMol EigenSolvers)
target_compile_definitions(MolDraw2D PRIVATE RDKIT_MOLDRAW2D_BUILD)

# if(RDK_BUILD_QT_DEMO && !RDK_BUILD_QT_SUPPORT) error
if (RDK_BUILD_QT_SUPPORT)
    add_subdirectory(Qt)
endif (RDK_BUILD_QT_SUPPORT)

if (RDK_BUILD_CAIRO_SUPPORT)
    find_package(Cairo REQUIRED)
    target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
    target_link_libraries(MolDraw2D PUBLIC Cairo::Cairo)
    target_sources(MolDraw2D PRIVATE MolDraw2DCairo.cpp DrawTextCairo.cpp)
    if (RDK_INSTALL_STATIC_LIBS AND NOT WIN32)
        target_compile_definitions(MolDraw2D_static PUBLIC "-DRDK_BUILD_CAIRO_SUPPORT")
        target_link_libraries(MolDraw2D_static PUBLIC Cairo::Cairo)
        target_sources(MolDraw2D_static PRIVATE MolDraw2DCairo.cpp DrawTextCairo.cpp)
    endif ()
    rdkit_headers(MolDraw2DCairo.h DEST GraphMol/MolDraw2D)
endif (RDK_BUILD_CAIRO_SUPPORT)

if (RDK_BUILD_FREETYPE_SUPPORT)
    target_compile_definitions(MolDraw2D PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
    target_sources(MolDraw2D PRIVATE DrawText.cpp
            DrawTextFT.cpp DrawTextFTSVG.cpp DrawTextFTJS.cpp)
    if (RDK_INSTALL_STATIC_LIBS AND NOT WIN32)
        target_compile_definitions(MolDraw2D_static PUBLIC "-DRDK_BUILD_FREETYPE_SUPPORT")
        target_sources(MolDraw2D_static PRIVATE DrawText.cpp
                DrawTextFT.cpp DrawTextFTSVG.cpp DrawTextFTJS.cpp)
    endif ()
    if (${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
        set(USE_FLAGS "-s USE_FREETYPE=1")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${USE_FLAGS}")
    endif ()
    find_package(Freetype REQUIRED)
    target_include_directories(MolDraw2D PUBLIC ${FREETYPE_INCLUDE_DIRS})
    target_link_libraries(MolDraw2D PUBLIC ${FREETYPE_LIBRARIES})
    if (RDK_INSTALL_STATIC_LIBS AND NOT WIN32)
        target_include_directories(MolDraw2D_static PUBLIC ${FREETYPE_INCLUDE_DIRS})
        target_link_libraries(MolDraw2D_static PUBLIC ${FREETYPE_LIBRARIES})
    endif ()
    if (RDK_BUILD_CAIRO_SUPPORT)
        target_sources(MolDraw2D PRIVATE DrawTextFTCairo.cpp)
        if (RDK_INSTALL_STATIC_LIBS AND NOT WIN32)
            target_sources(MolDraw2D_static PRIVATE DrawTextFTCairo.cpp)
        endif ()
    endif ()
endif (RDK_BUILD_FREETYPE_SUPPORT)

rdkit_test(moldraw2DTest1 test1.cpp LINK_LIBRARIES
        MolDraw2D)

rdkit_catch_test(moldraw2DTestCatch catch_tests.cpp LINK_LIBRARIES
        MolDraw2D CIPLabeler)

rdkit_test(moldraw2DRxnTest1 rxn_test1.cpp LINK_LIBRARIES
        MolDraw2D)

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
