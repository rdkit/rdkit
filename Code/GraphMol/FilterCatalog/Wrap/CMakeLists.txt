if(RDK_USE_BOOST_SERIALIZATION)
    set(T_LIBS ${Boost_LIBRARIES})
    find_package(Boost 1.39.0 COMPONENTS serialization)
    if (Boost_SERIALIZATION_LIBRARY)
      set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
    else()
      message("== Making Python FilterCatalog without boost Serialization support (MISSING)")
    endif()
else()
      message("== Making Python FilterCatalog without boost Serialization support")
endif()

rdkit_python_extension(rdfiltercatalog
                       rdfiltercatalog.cpp
                       FilterCatalog.cpp
                       DEST Chem 
                       LINK_LIBRARIES
                       FilterCatalog 
Subgraphs SubstructMatch SmilesParse Catalogs FileParsers GraphMol DataStructs
RDGeometryLib RDGeneral RDBoost  
${Boost_SERIALIZATION_LIBRARY} 
) 

add_pytest(pyFilterCatalog ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
