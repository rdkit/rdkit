project(catch_builder CXX)

# Includes Catch in the project:
include(ExternalProject)

ExternalProject_Add(
    catch
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch
    GIT_REPOSITORY https://github.com/philsquared/Catch.git
    GIT_TAG "v2.1.1"
    TIMEOUT 10
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch source_dir)
set(CATCH_MODULE_PATH ${source_dir} CACHE INTERNAL "location of Catch")
set(CATCH_INCLUDE_DIR ${CATCH_MODULE_PATH}/single_include CACHE INTERNAL "Path to include folder for Catch")
