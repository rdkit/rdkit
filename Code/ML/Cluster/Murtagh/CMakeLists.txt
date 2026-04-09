

rdkit_library(hc hc.c hcdriver.c)
target_compile_definitions(hc PRIVATE RDKIT_HC_BUILD)
if(RDK_USE_PRECOMPILED_HEADERS)
set_target_properties(hc PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)
if(RDK_INSTALL_STATIC_LIBS)
set_target_properties(hc_static PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)
endif(RDK_INSTALL_STATIC_LIBS)
endif(RDK_USE_PRECOMPILED_HEADERS)


set_target_properties(hc PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

rdkit_python_extension(Clustering Clustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral)



