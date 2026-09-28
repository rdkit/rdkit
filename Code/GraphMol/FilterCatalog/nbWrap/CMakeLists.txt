remove_definitions(-DRDKIT_FILTERCATALOG_BUILD)
rdkit_nanobind_extension(rdfiltercatalog
                         rdfiltercatalog.cpp
                         FilterCatalog.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         FilterCatalog)

add_pytest(pyFilterCatalog ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
