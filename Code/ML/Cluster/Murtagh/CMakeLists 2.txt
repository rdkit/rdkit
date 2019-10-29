
add_definitions(-DRDKIT_HC_BUILD)

rdkit_library(hc hc.c hcdriver.c)
set_target_properties(hc PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

rdkit_python_extension(Clustering Clustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral RDBoost)



