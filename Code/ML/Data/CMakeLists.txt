if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
rdkit_python_extension(cQuantize cQuantize.cpp
                       DEST ML/Data
                       LINK_LIBRARIES InfoTheory RDGeneral)
endif()

if(RDK_BUILD_NANOBIND_WRAPPERS)
  rdkit_nanobind_extension(cQuantize nbQuantize.cpp
                           DEST ML/Data
                           LINK_LIBRARIES InfoTheory RDGeneral)
endif()
