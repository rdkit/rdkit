if(RDK_USE_BOOST_SERIALIZATION AND Boost_SERIALIZATION_LIBRARY)
    set(RDKit_SERIALIZATION_LIBS ${Boost_SERIALIZATION_LIBRARY})
else()
    message("== Making FilterCatalog without boost Serialization support")
    set(RDKit_SERIALIZATION_LIBS )
endif()

rdkit_library(PortableSerializationLibs portable_binary_iarchive.cpp portable_binary_oarchive.cpp LINK_LIBRARIES RDGeneral ${RDKit_SERIALIZATION_LIBS})
target_compile_definitions(PortableSerializationLibs PRIVATE RDKIT_PORTABLEEIGENSERIALIZATITON_BUILD)
rdkit_library(EigenSerializer EigenSerializer.cpp LINK_LIBRARIES RDGeneral PortableSerializationLibs)
target_compile_definitions(EigenSerializer PRIVATE RDKIT_EIGENSERIALIZATITON_BUILD)

rdkit_headers(portable_binary_archive.hpp portable_binary_oarchive.hpp portable_binary_iarchive.hpp DEST Numerics/EigenSerializer)
rdkit_headers(EigenSerializer.h EigenBaseAddons.h DEST Numerics/EigenSerializer)
rdkit_catch_test(testEigenSerializer 
                 EigenSerializer_catch.cpp catch_main.cpp 
                 LINK_LIBRARIES RDGeneral PortableSerializationLibs EigenSerializer)
