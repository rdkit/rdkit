

# FIX: need code to grab the pubchem-align3d source

rdkit_library(pubchem_align3d ./pubchem-align3d/shape_functions1.cpp ./pubchem-align3d/shape_functions2.cpp ./pubchem-align3d/shape_neighbor.cpp)
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(pubchem_align3d PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()
rdkit_library(PubChemShape PubChemShape.cpp SHARED
  LINK_LIBRARIES pubchem_align3d SmilesParse SubstructMatch)
target_compile_definitions(PubChemShape PRIVATE RDKIT_PUBCHEMSHAPE_BUILD)

rdkit_headers(PubChemShape.h DEST GraphMol)

#add_executable(sdf_align sdf_align.cpp )
#target_link_libraries(sdf_align PubChemShape FileParsers)

rdkit_catch_test(shape_test test.cpp LINK_LIBRARIES PubChemShape FileParsers)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()




