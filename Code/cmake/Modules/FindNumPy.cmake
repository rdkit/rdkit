# adapted from http://voxel.jouy.inra.fr/darcs/contrib-itk/WrapITK/ExternalProjects/PyBuffer/FindNUMARRAY.cmake
# Try to find numarray python package
# Once done this will define
#
# PYTHON_NUMPY_FOUND        - system has numarray development package and it should be used
# PYTHON_NUMPY_INCLUDE_PATH  - directory where the arrayobject.h header file can be found
#
#

find_package(Python3 REQUIRED COMPONENTS Interpreter)
IF(Python3_EXECUTABLE)
  execute_process(
    COMMAND
    ${Python3_EXECUTABLE} -c "import numpy;print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE NUMPY_ERROR
  )
ENDIF(Python3_EXECUTABLE)

IF(NUMPY_ERROR)
  message(FATAL_ERROR "ERROR: The RDKit requires numpy. Numpy not found for Python ${${Python3_EXECUTABLE}}")
ENDIF()

FIND_PATH(PYTHON_NUMPY_INCLUDE_PATH numpy/arrayobject.h
  "${NUMPY_PATH}/"
  DOC "Directory where the arrayobject.h header file can be found. This file is part of the numarray package"
  )

IF(NOT PYTHON_NUMPY_INCLUDE_PATH)
  message(FATAL_ERROR "ERROR: The RDKit requires numpy. File numpy/arrayobject.h not found in ${NUMPY_PATH}")
ENDIF()

IF(PYTHON_NUMPY_INCLUDE_PATH)
  SET (PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy development package is available")
ENDIF(PYTHON_NUMPY_INCLUDE_PATH)
