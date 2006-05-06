# Microsoft Developer Studio Project File - Name="wrapA" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=wrapA - Win32 DebugPython
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "wrapA.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "wrapA.mak" CFG="wrapA - Win32 DebugPython"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "wrapA - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "wrapA - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "wrapA - Win32 DebugPython" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName "wrapA"
# PROP Scc_LocalPath "."
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "wrapA - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 1
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "$(RDBASE)/External/oelib" /I "$(RDBASE)/Code/GraphMol" /I "$(PYTHONHOME)\include" /I "$(RDBASE)/External/$(BOOSTBASE)" /I "$(RDBASE)/Code" /I "$(RDBASE)/External/vflib-2.0/include" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /FD /Zm400 /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /machine:I386
# ADD LINK32 boost_python.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib $(PYTHONVERS).lib BLAS.lib LAPACK.lib /nologo /subsystem:windows /dll /machine:I386 /nodefaultlib:"msvcprtdlib.lib" /out:"Release/moduleA.dll" /libpath:"$(PYTHONHOME)\libs" /libpath:"$(RDBASE)\External\Lapack\win32" /libpath:"$(RDBASE)/External/$(BOOSTBASE)/libs/python/build/bin/boost_python.dll/msvc/release/runtime-link-dynamic/"
# SUBTRACT LINK32 /verbose /pdb:none /incremental:yes
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=copy dll
PostBuild_Cmds=copy        Release\moduleA.dll      $(RDBASE)\Python\Chem
# End Special Build Tool

!ELSEIF  "$(CFG)" == "wrapA - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 1
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /GR /GX /Zi /Od /I "$(RDBASE)/External/oelib" /I "$(RDBASE)/Code/GraphMol" /I "$(PYTHONHOME)\include" /I "$(RDBASE)/External/$(BOOSTBASE)" /I "$(RDBASE)/Code" /I "$(RDBASE)/External/vflib-2.0/include" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /FD /GZ /Zm400 /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 boost_python_debug.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib $(PYTHONVERS).lib /nologo /subsystem:windows /dll /debug /machine:I386 /nodefaultlib:"msvcprtdlib.lib" /out:"moduleA.dll" /pdbtype:sept /libpath:"$(PYTHONHOME)\libs" /libpath:"$(RDBASE)/External/$(BOOSTBASE)/libs/python/build/bin-stage"
# SUBTRACT LINK32 /nodefaultlib

!ELSEIF  "$(CFG)" == "wrapA - Win32 DebugPython"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "wrapA_Win32_DebugPython"
# PROP BASE Intermediate_Dir "wrapA_Win32_DebugPython"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "wrapA_Win32_DebugPython"
# PROP Intermediate_Dir "wrapA_Win32_DebugPython"
# PROP Ignore_Export_Lib 1
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /GR /GX /Zi /Od /I "$(RDBASE)\external\boost_1_24_0" /I "$(RDBASE)/Code/GraphMol" /I "$(PYTHONHOME)\include" /I "$(RDBASE)/External/$(BOOSTBASE)" /I "$(RDBASE)/Code" /I "$(RDBASE)/External/vflib-2.0/include" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "TEST_EXPORTS" /D "BOOST_DEBUG_PYTHON" /FD /GZ /Zm400 /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 boost_python.lib, vflib.lib,libDatastructs.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib $(PYTHONVERS).lib BLAS.lib LAPACK.lib boost_python.lib /nologo /subsystem:windows /dll /debug /machine:I386 /nodefaultlib:"msvcprtdlib.lib" /out:"DebugPython/moduleA_d.dll" /pdbtype:sept /libpath:"c:\python20\activepython-2.0.0-202\src\core\pcbuild" /libpath:"$(BOOSTBASE)\libs\python\build\bin-stage" /libpath:"$(PYTHONHOME)\libs" /libpath:"$(RDBASE)\External\Lapack\win32" /libpath:"$(RDBASE)/External/$(BOOSTBASE)/libs/python/build/bin-stage"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=copy dll
PostBuild_Cmds=copy DebugPython\moduleA_d.dll $(RDBASE)\Python\Chem
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "wrapA - Win32 Release"
# Name "wrapA - Win32 Debug"
# Name "wrapA - Win32 DebugPython"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\wrapA.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
