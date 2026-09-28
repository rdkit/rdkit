rdkit_nanobind_extension(rdShapeAlign
                       cshapealign.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       PubChemShape )
add_pytest(pyShapeAlign
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/test_rdshapealign.py)
