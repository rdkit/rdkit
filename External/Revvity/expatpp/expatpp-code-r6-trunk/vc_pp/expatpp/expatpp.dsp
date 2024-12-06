# Microsoft Developer Studio Project File - Name="expatpp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=expatpp - Win32 ReleaseMT
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "expatpp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "expatpp.mak" CFG="expatpp - Win32 ReleaseMT"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "expatpp - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 DebugMTDLL" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 ReleaseMTDLL" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 DebugMT" (based on "Win32 (x86) Static Library")
!MESSAGE "expatpp - Win32 ReleaseMT" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "expatpp - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_Release"
# PROP BASE Intermediate_Dir "expatpp___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /I "\oofile\expatpp\expat\xmlparse" /I "\oofile\expatpp\expat\xmltok" /I "\oofile\expatpp\src_pp" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D XMLTOKAPI="" /D XMLPARSEAPI="" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /ML /W3 /GX /O2 /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 DebugMTDLL"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_DebugMTDLL"
# PROP BASE Intermediate_Dir "expatpp___Win32_DebugMTDLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugMTDLL"
# PROP Intermediate_Dir "DebugMTDLL"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MDd /W3 /GX /Z7 /Od /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 ReleaseMTDLL"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_ReleaseMTDLL"
# PROP BASE Intermediate_Dir "expatpp___Win32_ReleaseMTDLL"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "ReleaseMTDLL"
# PROP Intermediate_Dir "ReleaseMTDLL"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /ML /W3 /GX /O2 /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MD /W3 /GX /O2 /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 DebugMT"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_DebugMT"
# PROP BASE Intermediate_Dir "expatpp___Win32_DebugMT"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "DebugMT"
# PROP Intermediate_Dir "DebugMT"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MTd /W3 /GX /Z7 /Od /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "_DEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "expatpp - Win32 ReleaseMT"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "expatpp___Win32_ReleaseMT"
# PROP BASE Intermediate_Dir "expatpp___Win32_ReleaseMT"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "ReleaseMT"
# PROP Intermediate_Dir "ReleaseMT"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /ML /W3 /GX /O2 /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MT /W3 /GX /O2 /I "..\..\expat\lib\\" /D "XML_MIN_SIZE" /D "XML_WINLIB" /D "WIN32" /D "_WINDOWS" /D "NDEBUG" /D "COMPILED_FROM_DSP" /D "_LIB" /D "XML_STATIC" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09
# ADD RSC /l 0xc09
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "expatpp - Win32 Debug"
# Name "expatpp - Win32 Release"
# Name "expatpp - Win32 DebugMTDLL"
# Name "expatpp - Win32 ReleaseMTDLL"
# Name "expatpp - Win32 DebugMT"
# Name "expatpp - Win32 ReleaseMT"
# Begin Group "expat files"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\expat\lib\ascii.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\asciitab.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\expat.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\iasciitab.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\internal.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\latin1tab.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\macconfig.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\nametab.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\utf8tab.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\winconfig.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\xmlparse.c
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\xmlrole.c
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\xmlrole.h
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\xmltok.c
# End Source File
# Begin Source File

SOURCE=..\..\expat\lib\xmltok.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\src_pp\expatpp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\src_pp\expatpp.h
# End Source File
# End Target
# End Project
