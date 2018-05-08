remove_definitions(-DRDKIT_RDBOOST_BUILD)
rdkit_python_extension(rdBase
                       RDBase.cpp
                       LINK_LIBRARIES
                       RDGeneral RDBoost)

