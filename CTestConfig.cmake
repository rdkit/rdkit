# For some reason, MEMORYCHECK_SUPPRESSIONS_FILE is not being caught, so I hardcoded it here
SET(MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck --error-exitcode=13 --time-stamp=yes --num-callers=20 --gen-suppressions=all --leak-check=full --show-reachable=no --trace-children=yes --suppressions=${RDKit_SOURCE_DIR}/Code/cmake/rdkit_valgrind.suppressions")
