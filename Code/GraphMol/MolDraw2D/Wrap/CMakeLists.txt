remove_definitions(-DRDKIT_MOLDRAW2D_BUILD)

rdkit_python_extension(rdMolDraw2D
                       rdMolDraw2D.cpp
                       DEST Chem/Draw
                       LINK_LIBRARIES MolDraw2D )

# rdMolDraw2D.cpp #defines USE_BETTER_ENUMS and so is incompatible with the
# plain-enum precompiled header; compile it without the PCH (see BetterEnums.h).
if(RDK_USE_PRECOMPILED_HEADERS)
  set_source_files_properties(rdMolDraw2D.cpp PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
endif()

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2D.py)
