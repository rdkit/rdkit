remove_definitions(-DRDKIT_RDBOOST_BUILD)
rdkit_nanobind_extension(rdBase
                       RDBase.cpp
                       LINK_LIBRARIES
                       RDGeneral )

