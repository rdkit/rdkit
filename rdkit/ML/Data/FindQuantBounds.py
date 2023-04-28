#
#  Copyright (C) 2001  greg Landrum
#

from rdkit.Dbase import DbConnection
from rdkit.ML.Data import Quantize


def runIt(namesAndTypes, dbConnect, nBounds, resCol, typesToDo=['float']):
  results = map(lambda x: x[0], dbConnect.GetColumns(namesAndTypes[resCol][0]))
  nPossibleRes = max(results) + 1
  for cName, cType in namesAndTypes:
    if cType in typesToDo:
      dList = map(lambda x: x[0], dbConnect.GetColumns(cName))
      qDat = Quantize.FindVarMultQuantBounds(dList, nBounds, results, nPossibleRes)
      print(cName, qDat)


def Usage():
  import sys
  msg = """
  Usage: FindQuantBounds [-r res_col -n bounds_per_var -i] dbName tableName
  Optional Arguments:
    -r: specify the number of the result column
    -n: specify the number of bounds to attempt to find for each variable
    -i: also find vars for integer values
"""
  print(msg)
  sys.exit(-1)


if __name__ == '__main__':
  import getopt
  import sys

  try:
    args, extras = getopt.getopt(sys.argv[1:], 'n:r:i')
  except Exception:
    Usage()

  if len(extras) != 2:
    Usage()

  nBounds = 1
  typesToDo = ['float']
  includeInts = 0
  resCol = -1
  for arg, val in args:
    if arg == '-i':
      includeInts = 1
      typesToDo.append('integer')
    elif arg == '-n':
      try:
        nBounds = int(val)
      except ValueError:
        Usage()
    elif arg == '-r':
      try:
        resCol = int(val)
      except ValueError:
        Usage()

  dbName = extras[0]
  tableName = extras[1]
  dbConnect = DbConnection.DbConnect(dbName, tableName)
  namesAndTypes = dbConnect.GetColumnNamesAndTypes()
  runIt(namesAndTypes, dbConnect, nBounds, resCol, typesToDo)
