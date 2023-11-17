if(RDK_BUILD_PYTHON_STUBS)
  string(REGEX REPLACE "\\\\" "/" CMAKE_SOURCE_DIR_FWDSLASH ${CMAKE_SOURCE_DIR})
  string(REGEX REPLACE "\\\\" "/" CMAKE_CURRENT_BINARY_DIR_FWDSLASH ${CMAKE_CURRENT_BINARY_DIR})
  string(REGEX REPLACE "\\\\" "/" PYTHON_INSTDIR_FWDSLASH ${PYTHON_INSTDIR})
  string(CONCAT RUN_GEN_RDKIT_STUBS_PY "message("
    "\"-- Installing rdkit-stubs into the following directories: "
    "${CMAKE_SOURCE_DIR_FWDSLASH}, "
    "${PYTHON_INSTDIR_FWDSLASH}\")\n"
    "set (SEPARATOR \"=====================================================================\n\")\n"
    "set (COMMON_FILENAME ${CMAKE_CURRENT_BINARY_DIR_FWDSLASH}/gen_rdkit_stubs)\n"
    "execute_process("
    "COMMAND "
    "${PYTHON_EXECUTABLE} -m Scripts.gen_rdkit_stubs "
    "\"${CMAKE_SOURCE_DIR_FWDSLASH}\" \"${PYTHON_INSTDIR_FWDSLASH}\" "
    "WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR_FWDSLASH}\" "
    "OUTPUT_FILE \"\${COMMON_FILENAME}.out\" "
    "ERROR_FILE \"\${COMMON_FILENAME}.err\" "
    "RESULT_VARIABLE RETCODE)\n"
    "if (NOT \"\${RETCODE}\" EQUAL 0)\n"
    "  file(READ \"\${COMMON_FILENAME}.err\" ERROR_STRING)\n"
    "  message(\"\n\${SEPARATOR}Installation of rdkit-stubs failed.\n\${ERROR_STRING}\${SEPARATOR}\")\n"
    "endif()\n"
  )
  install(CODE "${RUN_GEN_RDKIT_STUBS_PY}" COMPONENT python)
endif(RDK_BUILD_PYTHON_STUBS)
