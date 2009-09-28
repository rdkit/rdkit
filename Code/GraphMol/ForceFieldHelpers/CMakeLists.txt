# moved library build directives to UFF (actual location of ForceFieldHelpers 
# modules.. might have been a bad idea).
add_subdirectory(UFF)

rdkit_test(testForceFieldHelpers testHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers FileParsers DistGeomHelpers EigenSolvers Alignment)

add_subdirectory(Wrap)

