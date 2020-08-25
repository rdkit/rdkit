

rdkit_library(hc hc.c hcdriver.c)
target_compile_definitions(hc PRIVATE RDKIT_HC_BUILD)

set_target_properties(hc PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

rdkit_python_extension(Clustering Clustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral)



