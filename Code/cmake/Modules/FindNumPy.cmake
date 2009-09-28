# adapted from http://voxel.jouy.inra.fr/darcs/contrib-itk/WrapITK/ExternalProjects/PyBuffer/FindNUMARRAY.cmake
# Try to find numarray python package
# Once done this will define
#
# PYTHON_NUMPY_FOUND        - system has numarray development package and it should be used
# PYTHON_NUMPY_INCLUDE_PATH  - directory where the arrayobject.h header file can be found
#
#

find_package(PythonInterp REQUIRED)
IF(PYTHON_EXECUTABLE)
  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/det_npp.py 
       "try: import numpy; print numpy.get_include()\nexcept: pass\n")
  EXEC_PROGRAM("${PYTHON_EXECUTABLE}"
    ARGS "\"${CMAKE_CURRENT_BINARY_DIR}/det_npp.py\""
    OUTPUT_VARIABLE NUMPY_PATH
  )
ENDIF(PYTHON_EXECUTABLE)

FIND_PATH(PYTHON_NUMPY_INCLUDE_PATH numpy/arrayobject.h
  "${NUMPY_PATH}/"
  DOC "Directory where the arrayobject.h header file can be found. This file is part of the numarray package"
  )

IF(PYTHON_NUMPY_INCLUDE_PATH)
  SET (PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy development package is available")
ENDIF(PYTHON_NUMPY_INCLUDE_PATH)
