find_package(SWIG 4.2 REQUIRED)
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
  "RascalMCES;GeneralizedSubstruct;Abbreviations;ScaffoldNetwork;MolHash;RGroupDecomposition;SubstructLibrary;TautomerQuery;"
  "MolEnumerator;"
  "MolStandardize;FilterCatalog;Catalogs;FMCS;MolDraw2D;FileParsers;SmilesParse;MarvinParser;"
  "Depictor;SubstructMatch;ChemReactions;Fingerprints;ChemTransforms;"
  "Subgraphs;GraphMol;DataStructs;Trajectory;Descriptors;"
  "PartialCharges;MolTransforms;DistGeomHelpers;DistGeometry;"
  "ForceFieldHelpers;ForceField;EigenSolvers;Optimizer;MolAlign;O3AAlign;"
  "Alignment;SimDivPickers;RDGeometryLib;RDStreams;RDGeneral;"
)
if(RDK_USE_URF)
  set(swigRDKitLibList "${swigRDKitLibList};${RDK_URF_LIBS};")
endif(RDK_USE_URF)
set(swigRDKitLibs "")
foreach(swigRDKitLib ${swigRDKitLibList})
  set(swigRDKitLibs "${swigRDKitLibs}${swigRDKitLib}${swigRDKitLibSuffix};")
endforeach()




set(RDKit_Wrapper_Libs ${swigRDKitLibs})

message("SUFFIX: ${swigRDKitLibSuffix}")
message("JAVA_LIBS: ${RDKit_Wrapper_Libs}")


if(RDK_BUILD_SWIG_JAVA_WRAPPER)
add_subdirectory(gmwrapper)
endif()

if(RDK_BUILD_SWIG_CSHARP_WRAPPER)
add_subdirectory(csharp_wrapper)
endif()
