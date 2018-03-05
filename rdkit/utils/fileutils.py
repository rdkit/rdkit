#
#  Copyright (C) 2000-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" utility functions to help work with files

"""


class NoMatchFoundError(RuntimeError):
  pass


def MoveToMatchingLine(inFile, matchStr, fullMatch=0):
  """ skip forward in a file until a given string is found

     **Arguments**

       - inFile: a file object (or anything supporting a _readline()_ method)

       - matchStr: the string to search for

       - fullMatch: if nonzero, _matchStr_ must match the entire line

     **Returns**

       the matching line

     **Notes:**

       - if _matchStr_ is not found in the file, a NoMatchFound exception
         will be raised

  """
  inLine = inFile.readline()
  matched = 0
  while not matched and inLine:
    idx = inLine.find(matchStr)
    if (fullMatch and idx == 0) or (not fullMatch and idx > -1):
      matched = 1
    else:
      inLine = inFile.readline()
  if matched:
    return inLine
  else:
    raise NoMatchFoundError(matchStr)
