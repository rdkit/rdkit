find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")


set(swigRDKitLibSuffix "")
if(RDK_SWIG_STATIC AND (NOT MSVC))
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
"FilterCatalog;Catalogs;FMCS;MolDraw2D;FileParsers;SmilesParse;"
"Depictor;SubstructMatch;ChemReactions;Fingerprints;ChemTransforms;"
"Subgraphs;GraphMol;DataStructs;Trajectory;Descriptors;"
"PartialCharges;MolTransforms;DistGeomHelpers;DistGeometry;"
"ForceFieldHelpers;ForceField;EigenSolvers;Optimizer;MolAlign;"
"Alignment;SimDivPickers;RDGeometryLib;RDGeneral;"
)
foreach(swigRDKitLib ${swigRDKitLibList})
  set(swigRDKitLibs "${swigRDKitLibs}${swigRDKitLib}${swigRDKitLibSuffix};")
endforeach()
set(swigRDKitLibs "${swigRDKitLibs}${Boost_SERIALIZATION_LIBRARY}")

set(RDKit_Wrapper_Libs ${swigRDKitLibs})

if(RDK_BUILD_SWIG_JAVA_WRAPPER)
add_subdirectory(gmwrapper)
endif()

if(RDK_BUILD_SWIG_CSHARP_WRAPPER)
add_subdirectory(csharp_wrapper)
endif()
