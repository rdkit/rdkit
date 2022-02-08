find_package(SWIG 3.0 EXACT REQUIRED)
include(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${RDKit_CodeDir})
include_directories(${Boost_INCLUDE_DIRS})

SET(CMAKE_SWIG_FLAGS "")


set(swigRDKitLibSuffix "")
if(RDK_SWIG_STATIC AND (NOT WIN32))
  set(swigRDKitLibSuffix "_static")
endif()
set(swigRDKitLibList "")
if(RDK_BUILD_AVALON_SUPPORT)
  set(swigRDKitLibList "${swigRDKitLibList}AvalonLib;avalon_clib;")
endif(RDK_BUILD_AVALON_SUPPORT)
if(RDK_BUILD_INCHI_SUPPORT)
  set(swigRDKitLibList "${swigRDKitLibList}RDInchiLib;${INCHI_LIBRARIES};")
endif(RDK_BUILD_INCHI_SUPPORT)
set(swigRDKitLibList "${swigRDKitLibList}"
  "Abbreviations;ScaffoldNetwork;MolHash;RGroupDecomposition;SubstructLibrary;TautomerQuery;"
  "MolEnumerator;"
  "MolStandardize;FilterCatalog;Catalogs;FMCS;MolDraw2D;FileParsers;SmilesParse;"
  "Depictor;SubstructMatch;ChemReactions;Fingerprints;ChemTransforms;"
  "Subgraphs;GraphMol;DataStructs;Trajectory;Descriptors;"
  "PartialCharges;MolTransforms;DistGeomHelpers;DistGeometry;"
  "ForceFieldHelpers;ForceField;EigenSolvers;Optimizer;MolAlign;O3AAlign;"
  "Alignment;SimDivPickers;RDGeometryLib;RDStreams;RDGeneral;"
)
#if(RDK_BUILD_COORDGEN_SUPPORT)
#  set(swigRDKitLibList "${swigRDKitLibList};${RDK_COORDGEN_LIBS}")
#endif(RDK_BUILD_COORDGEN_SUPPORT)
if(RDK_USE_URF)
  set(swigRDKitLibList "${swigRDKitLibList};${RDK_URF_LIBS};")
endif(RDK_USE_URF)
set(swigRDKitLibs "")
foreach(swigRDKitLib ${swigRDKitLibList})
  set(swigRDKitLibs "${swigRDKitLibs}${swigRDKitLib}${swigRDKitLibSuffix};")
endforeach()
set(swigRDKitLibs "${swigRDKitLibs}${Boost_SERIALIZATION_LIBRARY};")

if(RDK_USE_BOOST_IOSTREAMS)
ADD_DEFINITIONS("-DRDK_USE_BOOST_IOSTREAMS")
if(WIN32)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set(swigRDKitLibs "${swigRDKitLibs}${Boost_IOSTREAMS_LIBRARY};${Boost_SYSTEM_LIBRARY};")
  
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS zlib)
  if(Boost_ZLIB_FOUND)
    set(swigRDKitLibs "${swigRDKitLibs}${Boost_ZLIB_LIBRARY};")
  endif()
else()
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set(swigRDKitLibs "${swigRDKitLibs}${Boost_IOSTREAMS_LIBRARY};${Boost_SYSTEM_LIBRARY};")
endif()
endif(RDK_USE_BOOST_IOSTREAMS)



set(RDKit_Wrapper_Libs ${swigRDKitLibs})

message("SUFFIX: ${swigRDKitLibSuffix}")
message("JAVA_LIBS: ${RDKit_Wrapper_Libs}")


if(RDK_BUILD_SWIG_JAVA_WRAPPER)
add_subdirectory(gmwrapper)
endif()

if(RDK_BUILD_SWIG_CSHARP_WRAPPER)
add_subdirectory(csharp_wrapper)
endif()
