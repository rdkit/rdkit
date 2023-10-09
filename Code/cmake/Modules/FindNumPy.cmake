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
    ${Python3_EXECUTABLE} -c "try:\n import numpy\n print(numpy.get_include())\nexcept Exception:\n pass\n"
    OUTPUT_VARIABLE NUMPY_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
ENDIF(Python3_EXECUTABLE)

FIND_PATH(PYTHON_NUMPY_INCLUDE_PATH numpy/arrayobject.h
  "${NUMPY_PATH}/"
  DOC "Directory where the arrayobject.h header file can be found. This file is part of the numarray package"
  )

IF(PYTHON_NUMPY_INCLUDE_PATH)
  SET (PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy development package is available")
ENDIF(PYTHON_NUMPY_INCLUDE_PATH)
