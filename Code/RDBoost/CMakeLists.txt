rdkit_library(RDBoost Wrap.cpp
              LINK_LIBRARIES ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

rdkit_headers(Exceptions.h Wrap.h PySequenceHolder.h
              list_indexing_suite.hpp
              DEST RDBoost)

add_subdirectory(Wrap)
