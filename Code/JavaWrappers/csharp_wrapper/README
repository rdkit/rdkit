To build on Windows:
====================

Since cmake doesn't know anything about C#, there's an unfortunate
manual step involved in this. 

  - Make sure that the cmake configuration variable
    RDK_BUILD_SWIG_CSHARP_WRAPPER is set to ON. 
  - Run cmake to generate the solution file and open it in Visual
    Studio. 
  - Select the option to add an existing project and add
    $RDBASE/Code/JavaWrappers/csharp_wrapper/RDKit2DotNet.csproj 
  - Right click on the added project (named RDKit2DotNet) and add a
    dependency to RDKFuncs (this is the project that creates the C++
    dll that the C# project needs) 
  - Build the RDKit2DotNet project.

Your bin directory
($RDBASE/Code/JavaWrappers/csharp_wrapper/bin/Release if you did a
release build) now contains two DLLs: 
  - RDKFuncs.dll is the C++ dll containing the RDKit functionality 
  - RDKit2DotNet.dll contains the C# wrapper. 
To use the wrappers in your own projects, you should copy both dlls
into your project directory and add a reference to RDKit2DotNet.dll 

The directory RDKitCSharpTest contains a sample test project and some
code that makes very basic use of the wrapper functionality. 


To build on a linux system with mono installed:
==============================================

  - Make sure that the cmake configuration variable
    RDK_BUILD_SWIG_CSHARP_WRAPPER is set to ON. 
  - Run cmake and then do "make install"

This will create, in this directory, two files:
  - libRDKFuncs.{so,dynlib} (name depends on OS): the C++ dll
    containing the RDKit functionality
  - RDKit2DotNet.dll contains the C# wrapper.
  
Build and run the test code:
-------------------------------
gmcs -out:test.exe -addmodule:./RDKit2DotNet.dll test.cs
mono test.exe
