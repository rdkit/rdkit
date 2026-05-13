remove_definitions(-DRDKIT_SCAFFOLDNETWORK_BUILD)
rdkit_nanobind_extension(rdScaffoldNetwork
                         rdScaffoldNetwork.cpp
                         DEST Chem/Scaffolds
                         LINK_LIBRARIES ScaffoldNetwork
                                        GraphMol)

add_pytest(pyScaffoldNetwork ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testScaffoldNetwork.py)
if(RDK_USE_BOOST_SERIALIZATION)
  add_pytest(pyScaffoldNetworkPickling ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testPickleScaffoldNetwork.py)
endif()
