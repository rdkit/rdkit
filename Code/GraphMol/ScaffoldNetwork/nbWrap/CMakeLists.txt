remove_definitions(-DRDKIT_SCAFFOLDNETWORK_BUILD)
rdkit_python_extension(rdScaffoldNetwork
                       rdScaffoldNetwork.cpp
                       DEST Chem/Scaffolds 
                       LINK_LIBRARIES ScaffoldNetwork 
GraphMol ) 

add_pytest(pyScaffoldNetwork ${CMAKE_CURRENT_SOURCE_DIR}/testScaffoldNetwork.py)
if (RDK_USE_BOOST_SERIALIZATION)
    add_pytest(pyScaffoldNetworkPickling ${CMAKE_CURRENT_SOURCE_DIR}/testPickleScaffoldNetwork.py)
endif()
