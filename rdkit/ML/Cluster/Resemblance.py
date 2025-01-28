# $Id$
#
# Copyright (C) 2001-2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" code for dealing with resemblance (metric) matrices

    Here's how the matrices are stored:

     '[(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)...]  (row,col), col>row'

     or, alternatively the matrix can be drawn, with indices as:

       || - || 0 || 1 || 3
       || - || - || 2 || 4
       || - || - || - || 5
       || - || - || - || -

     the index of a given (row,col) pair is:
       '(col*(col-1))/2 + row'

"""

import numpy


def EuclideanDistance(inData):
  """returns the euclidean metricMat between the points in _inData_

    **Arguments**

     - inData: a Numeric array of data points

    **Returns**

       a Numeric array with the metric matrix.  See the module documentation
       for the format.


  """
  nObjs = len(inData)
  res = numpy.zeros((nObjs * (nObjs - 1) / 2), float)
  nSoFar = 0
  for col in range(1, nObjs):
    for row in range(col):
      t = inData[row] - inData[col]
      res[nSoFar] = sum(t * t)
      nSoFar += 1
  return numpy.sqrt(res)


def CalcMetricMatrix(inData, metricFunc):
  """ generates a metric matrix

    **Arguments**
     - inData is assumed to be a list of clusters (or anything with
       a GetPosition() method)

     - metricFunc is the function to be used to generate the matrix


    **Returns**

      the metric matrix as a Numeric array

  """
  #   nObjs = len(inData)
  #   res = []
  inData = map(lambda x: x.GetPosition(), inData)
  return metricFunc(inData)


def FindMinValInList(mat, nObjs, minIdx=None):
  """ finds the minimum value in a metricMatrix and returns it and its indices

    **Arguments**

     - mat: the metric matrix

     - nObjs: the number of objects to be considered

     - minIdx: the index of the minimum value (value, row and column still need
       to be calculated

    **Returns**

      a 3-tuple containing:

        1) the row
        2) the column
        3) the minimum value itself

    **Notes**

      -this probably ain't the speediest thing on earth

  """
  assert len(mat) == nObjs * (nObjs - 1) / 2, 'bad matrix length in FindMinValInList'
  if minIdx is None:
    minIdx = numpy.argmin(mat)

  nSoFar = 0
  col = 0
  while nSoFar <= minIdx:
    col = col + 1
    nSoFar += col

  row = minIdx - nSoFar + col
  return row, col, mat[minIdx]


def ShowMetricMat(metricMat, nObjs):
  """ displays a metric matrix

   **Arguments**

    - metricMat: the matrix to be displayed

    - nObjs: the number of objects to display

  """
  assert len(metricMat) == nObjs * (nObjs - 1) / 2, 'bad matrix length in FindMinValInList'
  for row in range(nObjs):
    for col in range(nObjs):
      if col <= row:
        print('   ---    ', end='')
      else:
        print('%10.6f' % metricMat[(col * (col - 1)) / 2 + row], end='')
    print()


methods = [
  ("Euclidean", EuclideanDistance, "Euclidean Distance"),
]

if __name__ == '__main__':
  m = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
  nObjs = 5
  for i in range(10):
    print(i, FindMinValInList(m, nObjs, i))
