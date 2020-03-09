# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" command line utility for growing composite models

**Usage**

  _GrowComposite [optional args] filename_

**Command Line Arguments**

  - -n *count*: number of new models to build

  - -C *pickle file name*:  name of file containing composite upon which to build.

  - --inNote *note*: note to be used in loading composite models from the database
      for growing

  - --balTable *table name*:  table from which to take the original data set
     (for balancing)

  - --balWeight *weight*: (between 0 and 1) weighting factor for the new data
     (for balancing). OR, *weight* can be a list of weights

  - --balCnt *count*: number of individual models in the balanced composite
     (for balancing)

  - --balH: use only the holdout set from the original data set in the balancing
     (for balancing)

  - --balT: use only the training set from the original data set in the balancing
     (for balancing)

  - -S: shuffle the original data set
     (for balancing)

  - -r: randomize the activities of the original data set
     (for balancing)

  - -N *note*: note to be attached to the grown composite when it's saved in the
     database

  - --outNote *note*: equivalent to -N

  - -o *filename*: name of an output file to hold the pickled composite after
     it has been grown.
     If multiple balance weights are used, the weights will be added to
     the filenames.

  - -L *limit*: provide an (integer) limit on individual model complexity

  - -d *database name*: instead of reading the data from a QDAT file,
     pull it from a database.  In this case, the _filename_ argument
     provides the name of the database table containing the data set.

  - -p *tablename*: store persistence data in the database
     in table *tablename*

  - -l: locks the random number generator to give consistent sets
     of training and hold-out data.  This is primarily intended
     for testing purposes.

  - -g: be less greedy when training the models.

  - -G *number*: force trees to be rooted at descriptor *number*.

  - -D: show a detailed breakdown of the composite model performance
     across the training and, when appropriate, hold-out sets.

  - -t *threshold value*: use high-confidence predictions for the final
     analysis of the hold-out data.

  - -q *list string*:  Add QuantTrees to the composite and use the list
     specified in *list string* as the number of target quantization
     bounds for each descriptor.  Don't forget to include 0's at the
     beginning and end of *list string* for the name and value fields.
     For example, if there are 4 descriptors and you want 2 quant bounds
     apiece, you would use _-q "[0,2,2,2,2,0]"_.
     Two special cases:
       1) If you would like to ignore a descriptor in the model building,
          use '-1' for its number of quant bounds.
       2) If you have integer valued data that should not be quantized
          further, enter 0 for that descriptor.

  - -V: print the version number and exit

"""


import sys
import time

import numpy

from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import CompositeRun
from rdkit.ML import ScreenComposite, BuildComposite
from rdkit.ML.Composite import AdjustComposite
from rdkit.ML.Data import DataUtils, SplitData
import pickle

_runDetails = CompositeRun.CompositeRun()

__VERSION_STRING = "0.5.0"

_verbose = 1


def message(msg):
  """ emits messages to _sys.stdout_
    override this in modules which import this one to redirect output

    **Arguments**

      - msg: the string to be displayed

  """
  if _verbose:
    sys.stdout.write('%s\n' % (msg))


def GrowIt(details, composite, progressCallback=None, saveIt=1, setDescNames=0, data=None):
  """ does the actual work of building a composite model

    **Arguments**

      - details:  a _CompositeRun.CompositeRun_ object containing details
        (options, parameters, etc.) about the run

      - composite: the composite model to grow

      - progressCallback: (optional) a function which is called with a single
        argument (the number of models built so far) after each model is built.

      - saveIt: (optional) if this is nonzero, the resulting model will be pickled
        and dumped to the filename specified in _details.outName_

      - setDescNames: (optional) if nonzero, the composite's _SetInputOrder()_ method
        will be called using the results of the data set's _GetVarNames()_ method;
        it is assumed that the details object has a _descNames attribute which
        is passed to the composites _SetDescriptorNames()_ method.  Otherwise
        (the default), _SetDescriptorNames()_ gets the results of _GetVarNames()_.

      - data: (optional) the data set to be used.  If this is not provided, the
        data set described in details will be used.

    **Returns**

      the enlarged composite model


  """
  details.rundate = time.asctime()

  if data is None:
    fName = details.tableName.strip()
    if details.outName == '':
      details.outName = fName + '.pkl'
    if details.dbName == '':
      data = DataUtils.BuildQuantDataSet(fName)
    elif details.qBounds != []:
      details.tableName = fName
      data = details.GetDataSet()
    else:
      data = DataUtils.DBToQuantData(  # Function no longer defined
        details.dbName, fName, quantName=details.qTableName, user=details.dbUser,
        password=details.dbPassword)

  seed = composite._randomSeed
  DataUtils.InitRandomNumbers(seed)
  if details.shuffleActivities == 1:
    DataUtils.RandomizeActivities(data, shuffle=1, runDetails=details)
  elif details.randomActivities == 1:
    DataUtils.RandomizeActivities(data, shuffle=0, runDetails=details)

  namedExamples = data.GetNamedData()
  trainExamples = namedExamples
  nExamples = len(trainExamples)
  message('Training with %d examples' % (nExamples))
  message('\t%d descriptors' % (len(trainExamples[0]) - 2))
  nVars = data.GetNVars()
  nPossibleVals = composite.nPossibleVals
  attrs = list(range(1, nVars + 1))

  if details.useTrees:
    from rdkit.ML.DecTree import CrossValidate, PruneTree
    if details.qBounds != []:
      from rdkit.ML.DecTree import BuildQuantTree
      builder = BuildQuantTree.QuantTreeBoot
    else:
      from rdkit.ML.DecTree import ID3
      builder = ID3.ID3Boot
    driver = CrossValidate.CrossValidationDriver
    pruner = PruneTree.PruneTree

    if setDescNames:
      composite.SetInputOrder(data.GetVarNames())
    composite.Grow(trainExamples, attrs, [0] + nPossibleVals, buildDriver=driver, pruner=pruner,
                   nTries=details.nModels, pruneIt=details.pruneIt, lessGreedy=details.lessGreedy,
                   needsQuantization=0, treeBuilder=builder, nQuantBounds=details.qBounds,
                   startAt=details.startAt, maxDepth=details.limitDepth,
                   progressCallback=progressCallback, silent=not _verbose)

  else:
    from rdkit.ML.Neural import CrossValidate
    driver = CrossValidate.CrossValidationDriver
    composite.Grow(trainExamples, attrs, [0] + nPossibleVals, nTries=details.nModels,
                   buildDriver=driver, needsQuantization=0)

  composite.AverageErrors()
  composite.SortModels()
  modelList, counts, avgErrs = composite.GetAllData()
  counts = numpy.array(counts)
  avgErrs = numpy.array(avgErrs)
  composite._varNames = data.GetVarNames()

  for i in range(len(modelList)):
    modelList[i].NameModel(composite._varNames)

  # do final statistics
  weightedErrs = counts * avgErrs
  averageErr = sum(weightedErrs) / sum(counts)
  devs = (avgErrs - averageErr)
  devs = devs * counts
  devs = numpy.sqrt(devs * devs)
  avgDev = sum(devs) / sum(counts)
  if _verbose:
    message('# Overall Average Error: %%% 5.2f, Average Deviation: %%% 6.2f' %
            (100. * averageErr, 100. * avgDev))

  if details.bayesModel:
    composite.Train(trainExamples, verbose=0)

  badExamples = []
  if not details.detailedRes:
    if _verbose:
      message('Testing all examples')
    wrong = BuildComposite.testall(composite, namedExamples, badExamples)
    if _verbose:
      message('%d examples (%% %5.2f) were misclassified' %
              (len(wrong), 100. * float(len(wrong)) / float(len(namedExamples))))
    _runDetails.overall_error = float(len(wrong)) / len(namedExamples)

  if details.detailedRes:
    if _verbose:
      message('\nEntire data set:')
    resTup = ScreenComposite.ShowVoteResults(
      range(data.GetNPts()), data, composite, nPossibleVals[-1], details.threshold)
    nGood, nBad, _, avgGood, avgBad, _, voteTab = resTup
    nPts = len(namedExamples)
    nClass = nGood + nBad
    _runDetails.overall_error = float(nBad) / nClass
    _runDetails.overall_correct_conf = avgGood
    _runDetails.overall_incorrect_conf = avgBad
    _runDetails.overall_result_matrix = repr(voteTab)
    nRej = nClass - nPts
    if nRej > 0:
      _runDetails.overall_fraction_dropped = float(nRej) / nPts

  return composite


def GetComposites(details):
  res = []
  if details.persistTblName and details.inNote:
    conn = DbConnect(details.dbName, details.persistTblName)
    mdls = conn.GetData(fields='MODEL', where="where note='%s'" % (details.inNote))
    for row in mdls:
      rawD = row[0]
      res.append(pickle.loads(str(rawD)))
  elif details.composFileName:
    res.append(pickle.load(open(details.composFileName, 'rb')))
  return res


def BalanceComposite(details, composite, data1=None, data2=None):
  """ balances the composite using the parameters provided in details

   **Arguments**

     - details a _CompositeRun.RunDetails_ object

     - composite: the composite model to be balanced

     - data1: (optional) if provided, this should be the
       data set used to construct the original models

     - data2: (optional) if provided, this should be the
       data set used to construct the new individual models

  """
  if not details.balCnt or details.balCnt > len(composite):
    return composite
  message("Balancing Composite")

  #
  # start by getting data set 1: which is the data set used to build the
  #  original models
  #
  if data1 is None:
    message("\tReading First Data Set")
    fName = details.balTable.strip()
    tmp = details.tableName
    details.tableName = fName
    dbName = details.dbName
    details.dbName = details.balDb
    data1 = details.GetDataSet()
    details.tableName = tmp
    details.dbName = dbName
  if data1 is None:
    return composite
  details.splitFrac = composite._splitFrac
  details.randomSeed = composite._randomSeed
  DataUtils.InitRandomNumbers(details.randomSeed)
  if details.shuffleActivities == 1:
    DataUtils.RandomizeActivities(data1, shuffle=1, runDetails=details)
  elif details.randomActivities == 1:
    DataUtils.RandomizeActivities(data1, shuffle=0, runDetails=details)
  namedExamples = data1.GetNamedData()
  if details.balDoHoldout or details.balDoTrain:
    trainIdx, testIdx = SplitData.SplitIndices(len(namedExamples), details.splitFrac, silent=1)
    trainExamples = [namedExamples[x] for x in trainIdx]
    testExamples = [namedExamples[x] for x in testIdx]
    if details.filterFrac != 0.0:
      trainIdx, temp = DataUtils.FilterData(trainExamples, details.filterVal, details.filterFrac,
                                            -1, indicesOnly=1)
      tmp = [trainExamples[x] for x in trainIdx]
      testExamples += [trainExamples[x] for x in temp]
      trainExamples = tmp
    if details.balDoHoldout:
      testExamples, trainExamples = trainExamples, testExamples
  else:
    trainExamples = namedExamples
  dataSet1 = trainExamples
  cols1 = [x.upper() for x in data1.GetVarNames()]
  data1 = None

  #
  # now grab data set 2: the data used to build the new individual models
  #
  if data2 is None:
    message("\tReading Second Data Set")
    data2 = details.GetDataSet()
  if data2 is None:
    return composite
  details.splitFrac = composite._splitFrac
  details.randomSeed = composite._randomSeed
  DataUtils.InitRandomNumbers(details.randomSeed)
  if details.shuffleActivities == 1:
    DataUtils.RandomizeActivities(data2, shuffle=1, runDetails=details)
  elif details.randomActivities == 1:
    DataUtils.RandomizeActivities(data2, shuffle=0, runDetails=details)
  dataSet2 = data2.GetNamedData()
  cols2 = [x.upper() for x in data2.GetVarNames()]
  data2 = None

  # and balance it:
  res = []
  weights = details.balWeight
  if not isinstance(weights, (tuple, list)):
    weights = (weights, )
  for weight in weights:
    message("\tBalancing with Weight: %.4f" % (weight))
    res.append(
      AdjustComposite.BalanceComposite(composite, dataSet1, dataSet2, weight, details.balCnt,
                                       names1=cols1, names2=cols2))
  return res


def ShowVersion(includeArgs=0):
  """ prints the version number

  """
  print('This is GrowComposite.py version %s' % (__VERSION_STRING))
  if includeArgs:
    print('command line was:')
    print(' '.join(sys.argv))


def Usage():
  """ provides a list of arguments for when this is used from the command line

  """
  print(__doc__)
  sys.exit(-1)


def SetDefaults(runDetails=None):
  """  initializes a details object with default values

      **Arguments**

        - details:  (optional) a _CompositeRun.CompositeRun_ object.
          If this is not provided, the global _runDetails will be used.

      **Returns**

        the initialized _CompositeRun_ object.


  """
  if runDetails is None:
    runDetails = _runDetails
  return CompositeRun.SetDefaults(runDetails)


def ParseArgs(runDetails):
  """ parses command line arguments and updates _runDetails_

      **Arguments**

        - runDetails:  a _CompositeRun.CompositeRun_ object.

  """
  import getopt
  args, extra = getopt.getopt(sys.argv[1:], 'P:o:n:p:b:sf:F:v:hlgd:rSTt:Q:q:DVG:L:C:N:',
                              ['inNote=',
                               'outNote=',
                               'balTable=',
                               'balWeight=',
                               'balCnt=',
                               'balH',
                               'balT',
                               'balDb=', ])
  runDetails.inNote = ''
  runDetails.composFileName = ''
  runDetails.balTable = ''
  runDetails.balWeight = (0.5, )
  runDetails.balCnt = 0
  runDetails.balDoHoldout = 0
  runDetails.balDoTrain = 0
  runDetails.balDb = ''
  for arg, val in args:
    if arg == '-n':
      runDetails.nModels = int(val)
    elif arg == '-C':
      runDetails.composFileName = val
    elif arg == '--balTable':
      runDetails.balTable = val
    elif arg == '--balWeight':
      runDetails.balWeight = eval(val)
      if not isinstance(runDetails.balWeight, (tuple, list)):
        runDetails.balWeight = (runDetails.balWeight, )
    elif arg == '--balCnt':
      runDetails.balCnt = int(val)
    elif arg == '--balH':
      runDetails.balDoHoldout = 1
    elif arg == '--balT':
      runDetails.balDoTrain = 1
    elif arg == '--balDb':
      runDetails.balDb = val
    elif arg == '--inNote':
      runDetails.inNote = val
    elif arg == '-N' or arg == '--outNote':
      runDetails.note = val
    elif arg == '-o':
      runDetails.outName = val
    elif arg == '-p':
      runDetails.persistTblName = val
    elif arg == '-r':
      runDetails.randomActivities = 1
    elif arg == '-S':
      runDetails.shuffleActivities = 1
    elif arg == '-h':
      Usage()
    elif arg == '-l':
      runDetails.lockRandom = 1
    elif arg == '-g':
      runDetails.lessGreedy = 1
    elif arg == '-G':
      runDetails.startAt = int(val)
    elif arg == '-d':
      runDetails.dbName = val
    elif arg == '-T':
      runDetails.useTrees = 0
    elif arg == '-t':
      runDetails.threshold = float(val)
    elif arg == '-D':
      runDetails.detailedRes = 1
    elif arg == '-L':
      runDetails.limitDepth = int(val)
    elif arg == '-q':
      qBounds = eval(val)
      assert isinstance(qBounds,
                        (tuple, list)), 'bad argument type for -q, specify a list as a string'
      runDetails.qBoundCount = val
      runDetails.qBounds = qBounds
    elif arg == '-Q':
      qBounds = eval(val)
      assert type(qBounds) in [type([]), type(
        ())], 'bad argument type for -Q, specify a list as a string'
      runDetails.activityBounds = qBounds
      runDetails.activityBoundsVals = val
    elif arg == '-V':
      ShowVersion()
      sys.exit(0)
    else:
      print('bad argument:', arg, file=sys.stderr)
      Usage()
  runDetails.tableName = extra[0]
  if not runDetails.balDb:
    runDetails.balDb = runDetails.dbName


if __name__ == '__main__':
  if len(sys.argv) < 2:
    Usage()

  _runDetails.cmd = ' '.join(sys.argv)
  SetDefaults(_runDetails)
  ParseArgs(_runDetails)

  ShowVersion(includeArgs=1)

  initModels = GetComposites(_runDetails)
  nModels = len(initModels)
  if nModels > 1:
    for i in range(nModels):
      sys.stderr.write(
        '---------------------------------\n\tDoing %d of %d\n---------------------------------\n' %
        (i + 1, nModels))
      composite = GrowIt(_runDetails, initModels[i], setDescNames=1)
      if _runDetails.balTable and _runDetails.balCnt:
        composites = BalanceComposite(_runDetails, composite)
      else:
        composites = [composite]
      for mdl in composites:
        mdl.ClearModelExamples()
      if _runDetails.outName:
        nWeights = len(_runDetails.balWeight)
        if nWeights == 1:
          outName = _runDetails.outName
          composites[0].Pickle(outName)
        else:
          for i in range(nWeights):
            weight = int(100 * _runDetails.balWeight[i])
            model = composites[i]
            outName = '%s.%d.pkl' % (_runDetails.outName.split('.pkl')[0], weight)
            model.Pickle(outName)
      if _runDetails.persistTblName and _runDetails.dbName:
        message('Updating results table %s:%s' % (_runDetails.dbName, _runDetails.persistTblName))
        if (len(_runDetails.balWeight)) > 1:
          message('WARNING: updating results table with models having different weights')
        # save the composite
        for i in range(len(composites)):
          _runDetails.model = pickle.dumps(composites[i])
          _runDetails.Store(db=_runDetails.dbName, table=_runDetails.persistTblName)
  elif nModels == 1:
    composite = GrowIt(_runDetails, initModels[0], setDescNames=1)
    if _runDetails.balTable and _runDetails.balCnt:
      composites = BalanceComposite(_runDetails, composite)
    else:
      composites = [composite]
    for mdl in composites:
      mdl.ClearModelExamples()
    if _runDetails.outName:
      nWeights = len(_runDetails.balWeight)
      if nWeights == 1:
        outName = _runDetails.outName
        composites[0].Pickle(outName)
      else:
        for i in range(nWeights):
          weight = int(100 * _runDetails.balWeight[i])
          model = composites[i]
          outName = '%s.%d.pkl' % (_runDetails.outName.split('.pkl')[0], weight)
          model.Pickle(outName)
    if _runDetails.persistTblName and _runDetails.dbName:
      message('Updating results table %s:%s' % (_runDetails.dbName, _runDetails.persistTblName))
      if (len(composites)) > 1:
        message('WARNING: updating results table with models having different weights')
      for i in range(len(composites)):
        _runDetails.model = pickle.dumps(composites[i])
        _runDetails.Store(db=_runDetails.dbName, table=_runDetails.persistTblName)
  else:
    message("No models found")
