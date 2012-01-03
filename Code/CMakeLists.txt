
add_subdirectory(RDGeneral)
if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(RDBoost)
endif(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(DataStructs)
add_subdirectory(Geometry)

add_subdirectory(Numerics)
add_subdirectory(ForceField)
add_subdirectory(DistGeom)
add_subdirectory(Catalogs)

add_subdirectory(GraphMol)
add_subdirectory(Query)

add_subdirectory(DataManip)
add_subdirectory(SimDivPickers)

add_subdirectory(ML)
add_subdirectory(ChemicalFeatures)

if(RDK_BUILD_SWIG_WRAPPERS)
  add_subdirectory(JavaWrappers)
endif()

