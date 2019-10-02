# $Id$
#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" command line utility to report on the contributions of descriptors to
tree-based composite models

Usage:  AnalyzeComposite [optional args] <models>

      <models>: file name(s) of pickled composite model(s)
        (this is the name of the db table if using a database)

    Optional Arguments:

      -n number: the number of levels of each model to consider

      -d dbname: the database from which to read the models

      -N Note: the note string to search for to pull models from the database

      -v: be verbose whilst screening
"""


import sys

import numpy

from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import ScreenComposite
from rdkit.ML.Data import Stats
from rdkit.ML.DecTree import TreeUtils, Tree
import pickle


__VERSION_STRING = "2.2.0"


def ProcessIt(composites, nToConsider=3, verbose=0):
  composite = composites[0]
  nComposites = len(composites)
  ns = composite.GetDescriptorNames()
  # nDesc = len(ns)-2
  if len(ns) > 2:
    globalRes = {}

    nDone = 1
    descNames = {}
    for composite in composites:
      if verbose > 0:
        print('#------------------------------------')
        print('Doing: ', nDone)
      nModels = len(composite)
      nDone += 1
      res = {}
      for i in range(len(composite)):
        model = composite.GetModel(i)
        if isinstance(model, Tree.TreeNode):
          levels = TreeUtils.CollectLabelLevels(model, {}, 0, nToConsider)
          TreeUtils.CollectDescriptorNames(model, descNames, 0, nToConsider)
          for descId in levels.keys():
            v = res.get(descId, numpy.zeros(nToConsider, numpy.float))
            v[levels[descId]] += 1. / nModels
            res[descId] = v
      for k in res:
        v = globalRes.get(k, numpy.zeros(nToConsider, numpy.float))
        v += res[k] / nComposites
        globalRes[k] = v
      if verbose > 0:
        for k in res.keys():
          name = descNames[k]
          strRes = ', '.join(['%4.2f' % x for x in res[k]])
          print('%s,%s,%5.4f' % (name, strRes, sum(res[k])))

        print()

    if verbose >= 0:
      print('# Average Descriptor Positions')
    retVal = []
    for k in globalRes:
      name = descNames[k]
      if verbose >= 0:
        strRes = ', '.join(['%4.2f' % x for x in globalRes[k]])
        print('%s,%s,%5.4f' % (name, strRes, sum(globalRes[k])))
      tmp = [name]
      tmp.extend(globalRes[k])
      tmp.append(sum(globalRes[k]))
      retVal.append(tmp)
    if verbose >= 0:
      print()
  else:
    retVal = []
  return retVal


def ErrorStats(conn, where, enrich=1):
  fields = ('overall_error,holdout_error,overall_result_matrix,' +
            'holdout_result_matrix,overall_correct_conf,overall_incorrect_conf,' +
            'holdout_correct_conf,holdout_incorrect_conf')
  try:
    data = conn.GetData(fields=fields, where=where)
  except Exception:
    import traceback
    traceback.print_exc()
    return None
  nPts = len(data)
  if not nPts:
    sys.stderr.write('no runs found\n')
    return None
  overall = numpy.zeros(nPts, numpy.float)
  overallEnrich = numpy.zeros(nPts, numpy.float)
  oCorConf = 0.0
  oInCorConf = 0.0
  holdout = numpy.zeros(nPts, numpy.float)
  holdoutEnrich = numpy.zeros(nPts, numpy.float)
  hCorConf = 0.0
  hInCorConf = 0.0
  overallMatrix = None
  holdoutMatrix = None
  for i in range(nPts):
    if data[i][0] is not None:
      overall[i] = data[i][0]
      oCorConf += data[i][4]
      oInCorConf += data[i][5]
    if data[i][1] is not None:
      holdout[i] = data[i][1]
      haveHoldout = 1
    else:
      haveHoldout = 0
    tmpOverall = 1. * eval(data[i][2])
    if enrich >= 0:
      overallEnrich[i] = ScreenComposite.CalcEnrichment(tmpOverall, tgt=enrich)
    if haveHoldout:
      tmpHoldout = 1. * eval(data[i][3])
      if enrich >= 0:
        holdoutEnrich[i] = ScreenComposite.CalcEnrichment(tmpHoldout, tgt=enrich)
    if overallMatrix is None:
      if data[i][2] is not None:
        overallMatrix = tmpOverall
      if haveHoldout and data[i][3] is not None:
        holdoutMatrix = tmpHoldout
    else:
      overallMatrix += tmpOverall
      if haveHoldout:
        holdoutMatrix += tmpHoldout
    if haveHoldout:
      hCorConf += data[i][6]
      hInCorConf += data[i][7]

  avgOverall = sum(overall) / nPts
  oCorConf /= nPts
  oInCorConf /= nPts
  overallMatrix /= nPts
  oSort = numpy.argsort(overall)
  oMin = overall[oSort[0]]
  overall -= avgOverall
  devOverall = numpy.sqrt(sum(overall**2) / (nPts - 1))
  res = {}
  res['oAvg'] = 100 * avgOverall
  res['oDev'] = 100 * devOverall
  res['oCorrectConf'] = 100 * oCorConf
  res['oIncorrectConf'] = 100 * oInCorConf
  res['oResultMat'] = overallMatrix
  res['oBestIdx'] = oSort[0]
  res['oBestErr'] = 100 * oMin

  if enrich >= 0:
    mean, dev = Stats.MeanAndDev(overallEnrich)
    res['oAvgEnrich'] = mean
    res['oDevEnrich'] = dev

  if haveHoldout:
    avgHoldout = sum(holdout) / nPts
    hCorConf /= nPts
    hInCorConf /= nPts
    holdoutMatrix /= nPts
    hSort = numpy.argsort(holdout)
    hMin = holdout[hSort[0]]
    holdout -= avgHoldout
    devHoldout = numpy.sqrt(sum(holdout**2) / (nPts - 1))
    res['hAvg'] = 100 * avgHoldout
    res['hDev'] = 100 * devHoldout
    res['hCorrectConf'] = 100 * hCorConf
    res['hIncorrectConf'] = 100 * hInCorConf
    res['hResultMat'] = holdoutMatrix
    res['hBestIdx'] = hSort[0]
    res['hBestErr'] = 100 * hMin
    if enrich >= 0:
      mean, dev = Stats.MeanAndDev(holdoutEnrich)
      res['hAvgEnrich'] = mean
      res['hDevEnrich'] = dev
  return res


def ShowStats(statD, enrich=1):
  statD = statD.copy()
  statD['oBestIdx'] = statD['oBestIdx'] + 1
  txt = """
# Error Statistics:
\tOverall: %(oAvg)6.3f%% (%(oDev)6.3f)  %(oCorrectConf)4.1f/%(oIncorrectConf)4.1f
\t\tBest: %(oBestIdx)d %(oBestErr)6.3f%%""" % (statD)
  if 'hAvg' in statD:
    statD['hBestIdx'] = statD['hBestIdx'] + 1
    txt += """
\tHoldout: %(hAvg)6.3f%% (%(hDev)6.3f)  %(hCorrectConf)4.1f/%(hIncorrectConf)4.1f
\t\tBest: %(hBestIdx)d %(hBestErr)6.3f%%
  """ % (statD)
  print(txt)
  print()
  print('# Results matrices:')
  print('\tOverall:')
  tmp = numpy.transpose(statD['oResultMat'])
  colCounts = sum(tmp)
  rowCounts = sum(tmp, 1)
  for i in range(len(tmp)):
    if rowCounts[i] == 0:
      rowCounts[i] = 1
    row = tmp[i]
    print('\t\t', end='')
    for j in range(len(row)):
      print('% 6.2f' % row[j], end='')
    print('\t| % 4.2f' % (100. * tmp[i, i] / rowCounts[i]))
  print('\t\t', end='')
  for i in range(len(tmp)):
    print('------', end='')
  print()
  print('\t\t', end='')
  for i in range(len(tmp)):
    if colCounts[i] == 0:
      colCounts[i] = 1
    print('% 6.2f' % (100. * tmp[i, i] / colCounts[i]), end='')
  print()
  if enrich > -1 and 'oAvgEnrich' in statD:
    print('\t\tEnrich(%d): %.3f (%.3f)' % (enrich, statD['oAvgEnrich'], statD['oDevEnrich']))

  if 'hResultMat' in statD:
    print('\tHoldout:')
    tmp = numpy.transpose(statD['hResultMat'])
    colCounts = sum(tmp)
    rowCounts = sum(tmp, 1)
    for i in range(len(tmp)):
      if rowCounts[i] == 0:
        rowCounts[i] = 1
      row = tmp[i]
      print('\t\t', end='')
      for j in range(len(row)):
        print('% 6.2f' % row[j], end='')
      print('\t| % 4.2f' % (100. * tmp[i, i] / rowCounts[i]))
    print('\t\t', end='')
    for i in range(len(tmp)):
      print('------', end='')
    print()
    print('\t\t', end='')
    for i in range(len(tmp)):
      if colCounts[i] == 0:
        colCounts[i] = 1
      print('% 6.2f' % (100. * tmp[i, i] / colCounts[i]), end='')
    print()
    if enrich > -1 and 'hAvgEnrich' in statD:
      print('\t\tEnrich(%d): %.3f (%.3f)' % (enrich, statD['hAvgEnrich'], statD['hDevEnrich']))

  return


def Usage():
  print(__doc__)
  sys.exit(-1)


if __name__ == "__main__":
  import getopt
  try:
    args, extras = getopt.getopt(sys.argv[1:], 'n:d:N:vX', ('skip',
                                                            'enrich=', ))
  except Exception:
    Usage()

  count = 3
  db = None
  note = ''
  verbose = 0
  skip = 0
  enrich = 1
  for arg, val in args:
    if arg == '-n':
      count = int(val) + 1
    elif arg == '-d':
      db = val
    elif arg == '-N':
      note = val
    elif arg == '-v':
      verbose = 1
    elif arg == '--skip':
      skip = 1
    elif arg == '--enrich':
      enrich = int(val)
  composites = []
  if db is None:
    for arg in extras:
      composite = pickle.load(open(arg, 'rb'))
      composites.append(composite)
  else:
    tbl = extras[0]
    conn = DbConnect(db, tbl)
    if note:
      where = "where note='%s'" % (note)
    else:
      where = ''
    if not skip:
      pkls = conn.GetData(fields='model', where=where)
      composites = []
      for pkl in pkls:
        pkl = str(pkl[0])
        comp = pickle.loads(pkl)
        composites.append(comp)

  if len(composites):
    ProcessIt(composites, count, verbose=verbose)
  elif not skip:
    print('ERROR: no composite models found')
    sys.exit(-1)

  if db:
    res = ErrorStats(conn, where, enrich=enrich)
    if res:
      ShowStats(res)
