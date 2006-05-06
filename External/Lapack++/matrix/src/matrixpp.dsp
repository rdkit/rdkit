# Microsoft Developer Studio Project File - Name="matrixpp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=matrixpp - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "matrixpp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "matrixpp.mak" CFG="matrixpp - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "matrixpp - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "matrixpp - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "matrixpp - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "../../include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "matrixpp - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "../../include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "matrixpp - Win32 Release"
# Name "matrixpp - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\bmd.cpp
# End Source File
# Begin Source File

SOURCE=.\gmc.cpp
# End Source File
# Begin Source File

SOURCE=.\gmd.cpp
# End Source File
# Begin Source File

SOURCE=.\gmf.cpp
# End Source File
# Begin Source File

SOURCE=.\gmi.cpp
# End Source File
# Begin Source File

SOURCE=.\gmli.cpp
# End Source File
# Begin Source File

SOURCE=.\index.cpp
# End Source File
# Begin Source File

SOURCE=.\lavc.cpp
# End Source File
# Begin Source File

SOURCE=.\lavd.cpp
# End Source File
# Begin Source File

SOURCE=.\lavi.cpp
# End Source File
# Begin Source File

SOURCE=.\lavli.cpp
# End Source File
# Begin Source File

SOURCE=.\ltgmd.cpp
# End Source File
# Begin Source File

SOURCE=.\spdbmd.cpp
# End Source File
# Begin Source File

SOURCE=.\spdmd.cpp
# End Source File
# Begin Source File

SOURCE=.\spdtrmd.cpp
# End Source File
# Begin Source File

SOURCE=.\sybmd.cpp
# End Source File
# Begin Source File

SOURCE=.\symd.cpp
# End Source File
# Begin Source File

SOURCE=.\sytrmd.cpp
# End Source File
# Begin Source File

SOURCE=.\trmd.cpp
# End Source File
# Begin Source File

SOURCE=.\ultgmd.cpp
# End Source File
# Begin Source File

SOURCE=.\utgmd.cpp
# End Source File
# Begin Source File

SOURCE=.\uutgmd.cpp
# End Source File
# Begin Source File

SOURCE=.\vc.cpp
# End Source File
# Begin Source File

SOURCE=.\vd.cpp
# End Source File
# Begin Source File

SOURCE=.\vf.cpp
# End Source File
# Begin Source File

SOURCE=.\vi.cpp
# End Source File
# Begin Source File

SOURCE=.\vli.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# End Target
# End Project
