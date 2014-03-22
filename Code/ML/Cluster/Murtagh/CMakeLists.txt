rdkit_library(hc hc.c hcdriver.c)

rdkit_python_extension(Clustering Clustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral RDBoost)



