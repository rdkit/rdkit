
# Search for the maeparser library or clone the sources from GitHub
macro(find_or_clone_maeparser)

    find_package(maeparser QUIET)

    if(maeparser_FOUND)

        message(STATUS "Found compiled maeparser library at ${maeparser_DIR}")

    elseif(NOT "${maeparser_DIR}" STREQUAL "")

        message(FATAL_ERROR "*** Failed to find a compiled instance of maeparser under "
                    "'${maeparser_DIR}'.")

    else()

        set(maeparser_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser-${MAEPARSER_VERSION}")

        if(NOT EXISTS "${maeparser_DIR}/maeparser/CMakeLists.txt")
            file(DOWNLOAD "https://github.com/schrodinger/maeparser/archive/${MAEPARSER_VERSION}.tar.gz"
                "${maeparser_DIR}/maeparser-${MAEPARSER_VERSION}.tar.gz")

            execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf "maeparser-${MAEPARSER_VERSION}.tar.gz"
                WORKING_DIRECTORY "${maeparser_DIR}")

            file(RENAME "${maeparser_DIR}/maeparser-${MAEPARSER_VERSION}" "${maeparser_DIR}/maeparser")
        endif()

        if(EXISTS "${maeparser_DIR}/maeparser/CMakeLists.txt")
            message(STATUS "Downloaded MaeParser '${MAEPARSER_VERSION}' to ${maeparser_DIR}.")
        else()
            message(FATAL_ERROR "Failed getting or unpacking Maeparser '${MAEPARSER_VERSION}'.")
        endif()

        add_subdirectory("${maeparser_DIR}/maeparser")

        set(maeparser_INCLUDE_DIRS "${maeparser_DIR}")
        set(maeparser_LIBRARIES maeparser)

    endif()

endmacro()