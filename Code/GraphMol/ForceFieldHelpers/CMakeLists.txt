# moved back library build directives here

rdkit_library(ForceFieldHelpers UFF/AtomTyper.cpp UFF/Builder.cpp
              MMFF/AtomTyper.cpp MMFF/Builder.cpp CrystalFF/TorsionAngleM6.cpp
              CrystalFF/TorsionPreferences.cpp
              LINK_LIBRARIES SmilesParse SubstructMatch ForceField)
target_compile_definitions(ForceFieldHelpers PRIVATE RDKIT_FORCEFIELDHELPERS_BUILD)

if(RDK_BUILD_ANI)
rdkit_library(ANIForceFieldHelpers ANI/Builder.cpp
              LINK_LIBRARIES SmilesParse SubstructMatch ForceField ANIForceField)
target_compile_definitions(ANIForceFieldHelpers PRIVATE ANIRDKIT_FORCEFIELDHELPERS_BUILD)
rdkit_headers(ANI/Builder.h ANI/ANI.h DEST GraphMol/ForceFieldHelpers/ANI)
endif(RDK_BUILD_ANI)
rdkit_headers(FFConvenience.h DEST GraphMol/ForceFieldHelpers)
rdkit_headers(UFF/AtomTyper.h
              UFF/Builder.h UFF/UFF.h DEST GraphMol/ForceFieldHelpers/UFF)
rdkit_headers(MMFF/AtomTyper.h
              MMFF/Builder.h MMFF/MMFF.h DEST GraphMol/ForceFieldHelpers/MMFF)
rdkit_headers(CrystalFF/TorsionAngleM6.h CrystalFF/TorsionPreferences.h
              DEST GraphMol/ForceFieldHelpers/CrystalFF)
rdkit_headers(ANI/Builder.h DEST GraphMol/ForceFieldHelpers/ANI)

add_subdirectory(MMFF)
add_subdirectory(UFF)
add_subdirectory(CrystalFF)
if(RDK_BUILD_ANI)
add_subdirectory(ANI)
endif(RDK_BUILD_ANI)
if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
