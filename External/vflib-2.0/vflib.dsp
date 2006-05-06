# Microsoft Developer Studio Project File - Name="vflib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=vflib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "vflib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "vflib.mak" CFG="vflib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "vflib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "vflib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "vflib - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /I "$(RDBASE)/External/vflib-2.0/include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "vflib - Win32 Debug"

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
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "$(RDBASE)/External/vflib-2.0/include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
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

# Name "vflib - Win32 Release"
# Name "vflib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\argedit.cpp
# End Source File
# Begin Source File

SOURCE=.\src\argloader.cpp
# End Source File
# Begin Source File

SOURCE=.\src\argraph.cpp
# End Source File
# Begin Source File

SOURCE=.\src\error.cpp
# End Source File
# Begin Source File

SOURCE=.\src\gene.cpp
# End Source File
# Begin Source File

SOURCE=.\src\gene_mesh.cpp
# End Source File
# Begin Source File

SOURCE=.\src\match.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sd_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sortnodes.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ull_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ull_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf2_mono_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf2_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf2_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf_mono_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\vf_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=.\src\xsubgraph.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\include\allocpool.h
# End Source File
# Begin Source File

SOURCE=.\include\argedit.h
# End Source File
# Begin Source File

SOURCE=.\include\argloader.h
# End Source File
# Begin Source File

SOURCE=.\include\argraph.h
# End Source File
# Begin Source File

SOURCE=.\include\dict.h
# End Source File
# Begin Source File

SOURCE=.\include\error.h
# End Source File
# Begin Source File

SOURCE=.\include\gene.h
# End Source File
# Begin Source File

SOURCE=.\include\gene_mesh.h
# End Source File
# Begin Source File

SOURCE=.\include\match.h
# End Source File
# Begin Source File

SOURCE=.\include\sd_state.h
# End Source File
# Begin Source File

SOURCE=.\src\sortnodes.h
# End Source File
# Begin Source File

SOURCE=.\include\state.h
# End Source File
# Begin Source File

SOURCE=.\include\ull_state.h
# End Source File
# Begin Source File

SOURCE=.\include\ull_sub_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf2_mono_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf2_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf2_sub_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf_mono_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf_state.h
# End Source File
# Begin Source File

SOURCE=.\include\vf_sub_state.h
# End Source File
# Begin Source File

SOURCE=.\include\xsubgraph.h
# End Source File
# End Group
# End Target
# End Project
