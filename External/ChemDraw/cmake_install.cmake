# Install script for directory: /Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/expatpp/expatpp-code-r6-trunk/expat/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.dylib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-Debug.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-MinSizeRel.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-RelWithDebInfo.cmake" OPTIONAL)
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitChemDraw.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitChemDraw.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitChemDraw.dylib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-Debug.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-MinSizeRel.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/ChemDraw.dir/install-cxx-module-bmi-RelWithDebInfo.cmake" OPTIONAL)
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitRDChemDrawLib.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug"
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitRDChemDrawLib.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release"
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitRDChemDrawLib.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel"
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
      "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitRDChemDrawLib.1.dylib"
      )
    foreach(file
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.2025.03.1pre.dylib"
        "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRDKitRDChemDrawLib.1.dylib"
        )
      if(EXISTS "${file}" AND
         NOT IS_SYMLINK "${file}")
        execute_process(COMMAND /Users/brian/miniforge3/bin/install_name_tool
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo"
          -delete_rpath "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_plac/lib"
          "${file}")
        if(CMAKE_INSTALL_DO_STRIP)
          execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
        endif()
      endif()
    endforeach()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Debug/libRDKitRDChemDrawLib.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/Release/libRDKitRDChemDrawLib.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/MinSizeRel/libRDKitRDChemDrawLib.dylib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/lib/RelWithDebInfo/libRDKitRDChemDrawLib.dylib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/RDChemDrawLib.dir/install-cxx-module-bmi-Debug.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/RDChemDrawLib.dir/install-cxx-module-bmi-Release.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/RDChemDrawLib.dir/install-cxx-module-bmi-MinSizeRel.cmake" OPTIONAL)
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    include("/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/CMakeFiles/RDChemDrawLib.dir/install-cxx-module-bmi-RelWithDebInfo.cmake" OPTIONAL)
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/rdkit/GraphMol" TYPE FILE FILES "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/chemdraw.h")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/brian/miniforge3/conda-bld/rdkit_1733161966177/work/External/Revvity/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
