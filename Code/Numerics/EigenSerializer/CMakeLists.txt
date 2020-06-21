if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making FilterCatalog without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

rdkit_library(EigenSerializer EigenSerializer.cpp LINK_LIBRARIES RDGeneral ${RDKit_SERIALIZATION_LIBS})

rdkit_headers(EigenSerializer.h EigenBaseAddons.h DEST Numerics/EigenSerializer)
rdkit_catch_test(testEigenSerializer 
                 EigenSerializer_catch.cpp catch_main.cpp 
                 LINK_LIBRARIES RDGeneral ${RDKit_SERIALIZATION_LIBS} EigenSerializer)