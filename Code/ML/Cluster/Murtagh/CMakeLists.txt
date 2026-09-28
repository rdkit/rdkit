

rdkit_library(hc hc.c hcdriver.c)
target_compile_definitions(hc PRIVATE RDKIT_HC_BUILD)

set_target_properties(hc PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)


if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
rdkit_python_extension(Clustering Clustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral)
endif(RDK_BUILD_BOOST_PYTHON_WRAPPERS)

if(RDK_BUILD_NANOBIND_WRAPPERS)
rdkit_nanobind_extension(Clustering nbClustering.cpp
                       DEST ML/Cluster
                       LINK_LIBRARIES
                       hc RDGeneral)
endif(RDK_BUILD_NANOBIND_WRAPPERS)

