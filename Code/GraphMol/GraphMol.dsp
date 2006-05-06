# Microsoft Developer Studio Project File - Name="GraphMol" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=GraphMol - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "GraphMol.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "GraphMol.mak" CFG="GraphMol - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "GraphMol - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "GraphMol - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "GraphMol - Win32 Release"

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
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "$(RDBASE)/External/$(BOOSTBASE)" /I "." /I "$(RDBASE)/Code" /I "$(RDBASE)\external\lapack++\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "GraphMol - Win32 Debug"

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
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "$(RDBASE)/External/$(BOOSTBASE)" /I "." /I "$(RDBASE)/Code" /I "$(RDBASE)\external\lapack++\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Ge /FD /GZ /c
# SUBTRACT CPP /Fr
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

# Name "GraphMol - Win32 Release"
# Name "GraphMol - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\AddHs.cpp
# End Source File
# Begin Source File

SOURCE=.\Aromaticity.cpp
# End Source File
# Begin Source File

SOURCE=.\Atom.cpp
# End Source File
# Begin Source File

SOURCE=.\atomic_data.cpp
# End Source File
# Begin Source File

SOURCE=.\AtomIterators.cpp
# End Source File
# Begin Source File

SOURCE=.\Bond.cpp
# End Source File
# Begin Source File

SOURCE=.\BondIterators.cpp
# End Source File
# Begin Source File

SOURCE=.\Canon.cpp
# End Source File
# Begin Source File

SOURCE=.\ConjugHybrid.cpp
# End Source File
# Begin Source File

SOURCE=.\FindRings.cpp
# End Source File
# Begin Source File

SOURCE=.\GraphMol.cpp
# End Source File
# Begin Source File

SOURCE=.\Invariant\Invariant.cpp
# End Source File
# Begin Source File

SOURCE=.\Kekulize.cpp
# End Source File
# Begin Source File

SOURCE=.\MolDiscriminators.cpp
# End Source File
# Begin Source File

SOURCE=.\MolOps.cpp
# End Source File
# Begin Source File

SOURCE=.\MolPickler.cpp
# End Source File
# Begin Source File

SOURCE=.\PeriodicTable.cpp
# End Source File
# Begin Source File

SOURCE=.\QueryAtom.cpp
# End Source File
# Begin Source File

SOURCE=.\QueryBond.cpp
# End Source File
# Begin Source File

SOURCE=.\QueryOps.cpp
# End Source File
# Begin Source File

SOURCE=.\RankAtoms.cpp
# End Source File
# Begin Source File

SOURCE=.\ROMol.cpp
# End Source File
# Begin Source File

SOURCE=.\RWMol.cpp
# End Source File
# Begin Source File

SOURCE=.\General\types.cpp
# End Source File
# Begin Source File

SOURCE=.\General\utils.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\Atom.h
# End Source File
# Begin Source File

SOURCE=.\atomic_data.h
# End Source File
# Begin Source File

SOURCE=.\AtomProps.h
# End Source File
# Begin Source File

SOURCE=.\Bond.h
# End Source File
# Begin Source File

SOURCE=.\BondProps.h
# End Source File
# Begin Source File

SOURCE=.\GraphMol.h
# End Source File
# Begin Source File

SOURCE=.\MolOps.h
# End Source File
# Begin Source File

SOURCE=.\PeriodicTable.h
# End Source File
# End Group
# End Target
# End Project
