remove_definitions(-DRDKIT_SCAFFOLDNETWORK_BUILD)
rdkit_nanobind_extension(rdScaffoldNetwork
                         rdScaffoldNetwork.cpp
                         DEST Chem/Scaffolds
                         LINK_LIBRARIES ScaffoldNetwork
                                        GraphMol)

add_pytest(pyScaffoldNetwork ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testScaffoldNetwork.py)
add_pytest(pyScaffoldNetworkPickling ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testPickleScaffoldNetwork.py)
add_pytest(pyNbScaffoldNetworkPickling ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testNbPickleScaffoldNetwork.py)
