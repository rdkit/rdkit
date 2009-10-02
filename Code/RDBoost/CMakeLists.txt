rdkit_library(RDBoost Wrap.cpp 
              LINK_LIBRARIES ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

add_subdirectory(Wrap)
