include(ExternalProject)
set(EIGEN3_VERSION "3.3.9")
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
ExternalProject_Add(Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG ${EIGEN3_VERSION}
  GIT_SHALLOW ON
  SOURCE_DIR ${EIGEN3_INCLUDE_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
set_target_properties(Eigen3::Eigen PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})
add_dependencies(Eigen3::Eigen Eigen)
