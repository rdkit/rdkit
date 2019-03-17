find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${RDKit_CodeDir})
include_directories(${Boost_INCLUDE_DIRS})

SET(CMAKE_SWIG_FLAGS "")


set(swigRDKitLibSuffix "")
if(RDK_SWIG_STATIC AND ((NOT MSVC) OR (MSVC AND RDK_INSTALL_DLLS_MSVC)))
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
  "RGroupDecomposition;SubstructLibrary;"
  "MolStandardize;FilterCatalog;Catalogs;FMCS;MolDraw2D;FileParsers;SmilesParse;"
  "Depictor;SubstructMatch;ChemReactions;Fingerprints;ChemTransforms;"
  "Subgraphs;GraphMol;DataStructs;Trajectory;Descriptors;"
  "PartialCharges;MolTransforms;DistGeomHelpers;DistGeometry;"
  "ForceFieldHelpers;ForceField;EigenSolvers;Optimizer;MolAlign;"
  "Alignment;SimDivPickers;RDGeometryLib;RDStreams;RDGeneral;"
)
#if(RDK_BUILD_COORDGEN_SUPPORT)
#  set(swigRDKitLibList "${swigRDKitLibList};${RDK_COORDGEN_LIBS}")
#endif(RDK_BUILD_COORDGEN_SUPPORT)
set(swigRDKitLibs "")
foreach(swigRDKitLib ${swigRDKitLibList})
  set(swigRDKitLibs "${swigRDKitLibs}${swigRDKitLib}${swigRDKitLibSuffix};")
endforeach()
set(swigRDKitLibs "${swigRDKitLibs}${Boost_SERIALIZATION_LIBRARY};")


if(WIN32)
  find_package(Boost 1.56.0 COMPONENTS system iostreams zlib REQUIRED)
  set(swigRDKitLibs "${swigRDKitLibs}${Boost_IOSTREAMS_LIBRARY};${Boost_ZLIB_LIBRARY};${Boost_SYSTEM_LIBRARY};")
else()
  find_package(Boost 1.56.0 COMPONENTS system iostreams REQUIRED)
  set(swigRDKitLibs "${swigRDKitLibs}${Boost_IOSTREAMS_LIBRARY};${Boost_SYSTEM_LIBRARY};")
endif()



set(RDKit_Wrapper_Libs ${swigRDKitLibs})

message("SUFFIX: ${swigRDKitLibSuffix}")
message("JAVA_LIBS: ${RDKit_Wrapper_Libs}")


if(RDK_BUILD_SWIG_JAVA_WRAPPER)
add_subdirectory(gmwrapper)
endif()

if(RDK_BUILD_SWIG_CSHARP_WRAPPER)
add_subdirectory(csharp_wrapper)
endif()
