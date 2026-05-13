remove_definitions(-DRDKIT_MOLSTANDARDIZE_BUILD)
rdkit_python_extension(rdMolStandardize rdMolStandardize.cpp Validate.cpp
	Charge.cpp Fragment.cpp Normalize.cpp Metal.cpp Tautomer.cpp Pipeline.cpp
                       DEST Chem/MolStandardize
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize )

add_pytest(pyMolStandardize ${CMAKE_CURRENT_SOURCE_DIR}/testMolStandardize.py)
