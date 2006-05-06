# Microsoft Developer Studio Project File - Name="SmilesParse" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=SmilesParse - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "SmilesParse.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "SmilesParse.mak" CFG="SmilesParse - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "SmilesParse - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "SmilesParse - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "SmilesParse - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_LIB" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "$(RDBASE)\code\graphmol" /I "$(RDBASE)\Code" /I "." /I "$(RDBASE)/External/$(BOOSTBASE)" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D YYDEBUG=1 /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "SmilesParse - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GR /GX /ZI /Od /I "$(RDBASE)/code/graphmol" /I "$(RDBASE)\Code" /I "." /I "$(RDBASE)/External/$(BOOSTBASE)" /D "WIN32" /D "_DEBUG" /D "_LIB" /D "_MBCS" /D YYDEBUG=1 /Ge /FD /GZ /c
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

# Name "SmilesParse - Win32 Release"
# Name "SmilesParse - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\InputFiller.cpp
# End Source File
# Begin Source File

SOURCE=.\lex.yysmarts.cpp
# End Source File
# Begin Source File

SOURCE=.\lex.yysmiles.cpp
# End Source File
# Begin Source File

SOURCE=.\smarts.tab.cpp
# End Source File
# Begin Source File

SOURCE=.\smiles.tab.cpp
# End Source File
# Begin Source File

SOURCE=.\SmilesParse.cpp
# End Source File
# Begin Source File

SOURCE=.\SmilesParseOps.cpp
# End Source File
# Begin Source File

SOURCE=.\SmilesWrite.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\InputFiller.h
# End Source File
# Begin Source File

SOURCE=.\smiles.tab.h
# PROP Ignore_Default_Tool 1
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\smarts.flex

!IF  "$(CFG)" == "SmilesParse - Win32 Release"

# Begin Custom Build
InputPath=.\smarts.flex

"lex.yysmarts.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/flex -Pyysmarts_ -olex.yy.tmp.cpp smarts.flex 
	c:/cygwin/bin/grep -v unistd lex.yy.tmp.cpp > lex.yysmarts.cpp 
	
# End Custom Build

!ELSEIF  "$(CFG)" == "SmilesParse - Win32 Debug"

# Begin Custom Build
InputPath=.\smarts.flex

"lex.yysmarts.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/flex -Pyysmarts_ -olex.yy.tmp.cpp smarts.flex 
	c:/cygwin/bin/grep -v unistd lex.yy.tmp.cpp > lex.yysmarts.cpp 
	
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\smarts.y

!IF  "$(CFG)" == "SmilesParse - Win32 Release"

# Begin Custom Build
InputPath=.\smarts.y

"smarts.tab.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/bison -p yysmarts_ -t -d -o smarts.tab.cpp smarts.y

# End Custom Build

!ELSEIF  "$(CFG)" == "SmilesParse - Win32 Debug"

# Begin Custom Build
InputPath=.\smarts.y

"smarts.tab.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/bison -p yysmarts_ -t -d -o smarts.tab.cpp smarts.y

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\smiles.flex

!IF  "$(CFG)" == "SmilesParse - Win32 Release"

# Begin Custom Build
InputPath=.\smiles.flex

"lex.yysmiles.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/flex -Pyysmiles_ -olex.yy.tmp.cpp smiles.flex 
	c:/cygwin/bin/grep -v unistd lex.yy.tmp.cpp > lex.yysmiles.cpp 
	
# End Custom Build

!ELSEIF  "$(CFG)" == "SmilesParse - Win32 Debug"

# Begin Custom Build
InputPath=.\smiles.flex

"lex.yysmiles.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/flex -Pyysmiles_ -olex.yy.tmp.cpp smiles.flex 
	c:/cygwin/bin/grep -v unistd lex.yy.tmp.cpp > lex.yysmiles.cpp 
	
# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\smiles.y

!IF  "$(CFG)" == "SmilesParse - Win32 Release"

# Begin Custom Build
InputPath=.\smiles.y

"smiles.tab.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/bison -p yysmiles_ -t -d -o smiles.tab.cpp smiles.y

# End Custom Build

!ELSEIF  "$(CFG)" == "SmilesParse - Win32 Debug"

# Begin Custom Build
InputPath=.\smiles.y

"smiles.tab.cpp" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	c:/cygwin/bin/bison -p yysmiles_ -t -d -o smiles.tab.cpp smiles.y

# End Custom Build

!ENDIF 

# End Source File
# End Target
# End Project
