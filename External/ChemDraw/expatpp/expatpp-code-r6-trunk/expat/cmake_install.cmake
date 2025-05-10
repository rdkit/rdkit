# Install script for directory: /Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

set(CMAKE_BINARY_DIR "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work")

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Debug/libexpat.dylib")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -id "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Release/libexpat.dylib")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -id "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/MinSizeRel/libexpat.dylib")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -id "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/RelWithDebInfo/libexpat.dylib")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -id "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libexpat.dylib")
      endif()
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/lib/expat.h"
    "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/lib/expat_external.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/expat.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Debug/xmlwf")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -change "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Debug/libexpat.dylib" "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Release/xmlwf")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -change "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/Release/libexpat.dylib" "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/MinSizeRel/xmlwf")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -change "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/MinSizeRel/libexpat.dylib" "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      endif()
    endif()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/RelWithDebInfo/xmlwf")
    if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf" AND
       NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      execute_process(COMMAND "/Users/brian/miniforge3/bin/install_name_tool"
        -change "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/RelWithDebInfo/libexpat.dylib" "libexpat.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/xmlwf")
      endif()
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/CMakeFiles/xmlwf.dir/install-cxx-module-bmi-Debug.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/CMakeFiles/xmlwf.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/CMakeFiles/xmlwf.dir/install-cxx-module-bmi-MinSizeRel.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/CMakeFiles/xmlwf.dir/install-cxx-module-bmi-RelWithDebInfo.cmake" OPTIONAL)
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1" TYPE FILE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/doc/xmlwf.1")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
