# $Id$
#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" command line utility for building composite models

#DOC

**Usage**

  BuildComposite [optional args] filename

Unless indicated otherwise (via command line arguments), _filename_ is
a QDAT file.

**Command Line Arguments**

  - -o *filename*: name of the output file for the pickled composite

  - -n *num*: number of separate models to add to the composite

  - -p *tablename*: store persistence data in the database
     in table *tablename*

  - -N *note*: attach some arbitrary text to the persistence data

  - -b *filename*: name of the text file to hold examples from the
     holdout set which are misclassified

  - -s: split the data into training and hold-out sets before building
     the composite

  - -f *frac*: the fraction of data to use in the training set when the
     data is split

  - -r: randomize the activities (for testing purposes).  This ignores
     the initial distribution of activity values and produces each
     possible activity value with equal likliehood.

  - -S: shuffle the activities (for testing purposes) This produces
     a permutation of the input activity values.

  - -l: locks the random number generator to give consistent sets
     of training and hold-out data.  This is primarily intended
     for testing purposes.

  - -B: use a so-called Bayesian composite model.

  - -d *database name*: instead of reading the data from a QDAT file,
     pull it from a database.  In this case, the _filename_ argument
     provides the name of the database table containing the data set.

  - -D: show a detailed breakdown of the composite model performance
     across the training and, when appropriate, hold-out sets.

  - -P *pickle file name*: write out the pickled data set to the file

  - -F *filter frac*: filters the data before training to change the
     distribution of activity values in the training set.  *filter
     frac* is the fraction of the training set that should have the
     target value.  **See note below on data filtering.**

  - -v *filter value*: filters the data before training to change the
     distribution of activity values in the training set. *filter
     value* is the target value to use in filtering.  **See note below
     on data filtering.**

  - --modelFiltFrac *model filter frac*: Similar to filter frac above,
     in this case the data is filtered for each model in the composite
     rather than a single overall filter for a composite. *model
     filter frac* is the fraction of the training set for each model
     that should have the target value (*model filter value*).

  - --modelFiltVal *model filter value*: target value to use for
     filtering data before training each model in the composite.

  - -t *threshold value*: use high-confidence predictions for the
     final analysis of the hold-out data.

  - -Q *list string*: the values of quantization bounds for the
     activity value.  See the _-q_ argument for the format of *list
     string*.

  - --nRuns *count*: build *count* composite models

  - --prune: prune any models built

  - -h: print a usage message and exit.

  - -V: print the version number and exit

  *-*-*-*-*-*-*-*- Tree-Related Options -*-*-*-*-*-*-*-*

  - -g: be less greedy when training the models.

  - -G *number*: force trees to be rooted at descriptor *number*.

  - -L *limit*: provide an (integer) limit on individual model
     complexity

  - -q *list string*: Add QuantTrees to the composite and use the list
     specified in *list string* as the number of target quantization
     bounds for each descriptor.  Don't forget to include 0's at the
     beginning and end of *list string* for the name and value fields.
     For example, if there are 4 descriptors and you want 2 quant
     bounds apiece, you would use _-q "[0,2,2,2,2,0]"_.
     Two special cases:
       1) If you would like to ignore a descriptor in the model
          building, use '-1' for its number of quant bounds.
       2) If you have integer valued data that should not be quantized
          further, enter 0 for that descriptor.

  - --recycle: allow descriptors to be used more than once in a tree

  - --randomDescriptors=val: toggles growing random forests with val
      randomly-selected descriptors available at each node.


  *-*-*-*-*-*-*-*- KNN-Related Options -*-*-*-*-*-*-*-*

  - --doKnn: use K-Nearest Neighbors models

  - --knnK=*value*: the value of K to use in the KNN models

  - --knnTanimoto: use the Tanimoto metric in KNN models

  - --knnEuclid: use a Euclidean metric in KNN models

  *-*-*-*-*-*-*- Naive Bayes Classifier Options -*-*-*-*-*-*-*-*
  - --doNaiveBayes : use Naive Bayes classifiers

  - --mEstimateVal : the value to be used in the m-estimate formula
      If this is greater than 0.0, we use it to compute the conditional
      probabilities by the m-estimate

  *-*-*-*-*-*-*-*- SVM-Related Options -*-*-*-*-*-*-*-*

  **** NOTE: THESE ARE DISABLED ****

# #   - --doSVM: use Support-vector machines

# #   - --svmKernel=*kernel*: choose the type of kernel to be used for
# #     the SVMs.  Options are:
# #     The default is:

# #   - --svmType=*type*: choose the type of support-vector machine
# #     to be used.  Options are:
# #     The default is:

# #   - --svmGamma=*gamma*: provide the gamma value for the SVMs.  If this
# #     is not provided, a grid search will be carried out to determine an
# #     optimal *gamma* value for each SVM.

# #   - --svmCost=*cost*: provide the cost value for the SVMs.  If this is
# #     not provided, a grid search will be carried out to determine an
# #     optimal *cost* value for each SVM.

# #   - --svmWeights=*weights*: provide the weight values for the
# #     activities.  If provided this should be a sequence of (label,
# #     weight) 2-tuples *nActs* long.  If not provided, a weight of 1
# #     will be used for each activity.

# #   - --svmEps=*epsilon*: provide the epsilon value used to determine
# #     when the SVM has converged.  Defaults to 0.001

# #   - --svmDegree=*degree*: provide the degree of the kernel (when
# #     sensible) Defaults to 3

# #   - --svmCoeff=*coeff*: provide the coefficient for the kernel (when
# #     sensible) Defaults to 0

# #   - --svmNu=*nu*: provide the nu value for the kernel (when sensible)
# #     Defaults to 0.5

# #   - --svmDataType=*float*: if the data is contains only 1 and 0 s, specify by
# #     using binary. Defaults to float

# #   - --svmCache=*cache*: provide the size of the memory cache (in MB)
# #     to be used while building the SVM.  Defaults to 40

**Notes**

  - *Data filtering*: When there is a large disparity between the
    numbers of points with various activity levels present in the
    training set it is sometimes desirable to train on a more
    homogeneous data set.  This can be accomplished using filtering.
    The filtering process works by selecting a particular target
    fraction and target value.  For example, in a case where 95% of
    the original training set has activity 0 and ony 5% activity 1, we
    could filter (by randomly removing points with activity 0) so that
    30% of the data set used to build the composite has activity 1.


"""


import sys
import time

import numpy

from rdkit import DataStructs
from rdkit.Dbase import DbModule
from rdkit.ML import CompositeRun
from rdkit.ML import ScreenComposite
from rdkit.ML.Composite import Composite, BayesComposite
from rdkit.ML.Data import DataUtils, SplitData
from rdkit.utils import listutils
import pickle

# # from ML.SVM import SVMClassificationModel as SVM
_runDetails = CompositeRun.CompositeRun()

__VERSION_STRING = "3.2.3"

_verbose = 1


def message(msg):
  """ emits messages to _sys.stdout_
    override this in modules which import this one to redirect output

    **Arguments**

      - msg: the string to be displayed

  """
  if _verbose:
    sys.stdout.write('%s\n' % (msg))


def testall(composite, examples, badExamples=[]):
  """ screens a number of examples past a composite

    **Arguments**

      - composite: a composite model

      - examples: a list of examples (with results) to be screened

      - badExamples: a list to which misclassified examples are appended

    **Returns**

      a list of 2-tuples containing:

        1) a vote

        2) a confidence

      these are the votes and confidence levels for **misclassified** examples

  """
  wrong = []
  for example in examples:
    if composite.GetActivityQuantBounds():
      answer = composite.QuantizeActivity(example)[-1]
    else:
      answer = example[-1]
    res, conf = composite.ClassifyExample(example)
    if res != answer:
      wrong.append((res, conf))
      badExamples.append(example)

  return wrong


def GetCommandLine(details):
  """ #DOC

  """
  args = ['BuildComposite']
  args.append('-n %d' % (details.nModels))
  if details.filterFrac != 0.0:
    args.append('-F %.3f -v %d' % (details.filterFrac, details.filterVal))
  if details.modelFilterFrac != 0.0:
    args.append('--modelFiltFrac=%.3f --modelFiltVal=%d' % (details.modelFilterFrac,
                                                            details.modelFilterVal))
  if details.splitRun:
    args.append('-s -f %.3f' % (details.splitFrac))
  if details.shuffleActivities:
    args.append('-S')
  if details.randomActivities:
    args.append('-r')
  if details.threshold > 0.0:
    args.append('-t %.3f' % (details.threshold))
  if details.activityBounds:
    args.append('-Q "%s"' % (details.activityBoundsVals))
  if details.dbName:
    args.append('-d %s' % (details.dbName))
  if details.detailedRes:
    args.append('-D')
  if hasattr(details, 'noScreen') and details.noScreen:
    args.append('--noScreen')
  if details.persistTblName and details.dbName:
    args.append('-p %s' % (details.persistTblName))
  if details.note:
    args.append('-N %s' % (details.note))
  if details.useTrees:
    if details.limitDepth > 0:
      args.append('-L %d' % (details.limitDepth))
    if details.lessGreedy:
      args.append('-g')
    if details.qBounds:
      shortBounds = listutils.CompactListRepr(details.qBounds)
      if details.qBounds:
        args.append('-q "%s"' % (shortBounds))
    else:
      if details.qBounds:
        args.append('-q "%s"' % (details.qBoundCount))

    if details.pruneIt:
      args.append('--prune')
    if details.startAt:
      args.append('-G %d' % details.startAt)
    if details.recycleVars:
      args.append('--recycle')
    if details.randomDescriptors:
      args.append('--randomDescriptors=%d' % details.randomDescriptors)
  if details.useSigTrees:
    args.append('--doSigTree')
    if details.limitDepth > 0:
      args.append('-L %d' % (details.limitDepth))
    if details.randomDescriptors:
      args.append('--randomDescriptors=%d' % details.randomDescriptors)

  if details.useKNN:
    args.append('--doKnn --knnK %d' % (details.knnNeighs))
    if details.knnDistFunc == 'Tanimoto':
      args.append('--knnTanimoto')
    else:
      args.append('--knnEuclid')

  if details.useNaiveBayes:
    args.append('--doNaiveBayes')
    if details.mEstimateVal >= 0.0:
      args.append('--mEstimateVal=%.3f' % details.mEstimateVal)

  # #   if details.useSVM:
  # #     args.append('--doSVM')
  # #     if details.svmKernel:
  # #       for k in SVM.kernels.keys():
  # #         if SVM.kernels[k]==details.svmKernel:
  # #           args.append('--svmKernel=%s'%k)
  # #           break
  # #     if details.svmType:
  # #       for k in SVM.machineTypes.keys():
  # #         if SVM.machineTypes[k]==details.svmType:
  # #           args.append('--svmType=%s'%k)
  # #           break
  # #     if details.svmGamma:
  # #       args.append('--svmGamma=%f'%details.svmGamma)
  # #     if details.svmCost:
  # #       args.append('--svmCost=%f'%details.svmCost)
  # #     if details.svmWeights:
  # #       args.append("--svmWeights='%s'"%str(details.svmWeights))
  # #     if details.svmDegree:
  # #       args.append('--svmDegree=%d'%details.svmDegree)
  # #     if details.svmCoeff:
  # #       args.append('--svmCoeff=%d'%details.svmCoeff)
  # #     if details.svmEps:
  # #       args.append('--svmEps=%f'%details.svmEps)
  # #     if details.svmNu:
  # #       args.append('--svmNu=%f'%details.svmNu)
  # #     if details.svmCache:
  # #       args.append('--svmCache=%d'%details.svmCache)
  # #     if detail.svmDataType:
  # #       args.append('--svmDataType=%s'%details.svmDataType)
  # #     if not details.svmShrink:
  # #       args.append('--svmShrink')

  if details.replacementSelection:
    args.append('--replacementSelection')

  # this should always be last:
  if details.tableName:
    args.append(details.tableName)

  return ' '.join(args)


def RunOnData(details, data, progressCallback=None, saveIt=1, setDescNames=0):
  if details.lockRandom:
    seed = details.randomSeed
  else:
    import random
    seed = (random.randint(0, 1e6), random.randint(0, 1e6))
  DataUtils.InitRandomNumbers(seed)
  testExamples = []
  if details.shuffleActivities == 1:
    DataUtils.RandomizeActivities(data, shuffle=1, runDetails=details)
  elif details.randomActivities == 1:
    DataUtils.RandomizeActivities(data, shuffle=0, runDetails=details)

  namedExamples = data.GetNamedData()
  if details.splitRun == 1:
    trainIdx, testIdx = SplitData.SplitIndices(
      len(namedExamples), details.splitFrac, silent=not _verbose)

    trainExamples = [namedExamples[x] for x in trainIdx]
    testExamples = [namedExamples[x] for x in testIdx]
  else:
    testExamples = []
    testIdx = []
    trainIdx = list(range(len(namedExamples)))
    trainExamples = namedExamples

  if details.filterFrac != 0.0:
    # if we're doing quantization on the fly, we need to handle that here:
    if hasattr(details, 'activityBounds') and details.activityBounds:
      tExamples = []
      bounds = details.activityBounds
      for pt in trainExamples:
        pt = pt[:]
        act = pt[-1]
        placed = 0
        bound = 0
        while not placed and bound < len(bounds):
          if act < bounds[bound]:
            pt[-1] = bound
            placed = 1
          else:
            bound += 1
        if not placed:
          pt[-1] = bound
        tExamples.append(pt)
    else:
      bounds = None
      tExamples = trainExamples
    trainIdx, temp = DataUtils.FilterData(tExamples, details.filterVal, details.filterFrac, -1,
                                          indicesOnly=1)
    tmp = [trainExamples[x] for x in trainIdx]
    testExamples += [trainExamples[x] for x in temp]
    trainExamples = tmp

    counts = DataUtils.CountResults(trainExamples, bounds=bounds)
    ks = counts.keys()
    ks.sort()
    message('Result Counts in training set:')
    for k in ks:
      message(str((k, counts[k])))
    counts = DataUtils.CountResults(testExamples, bounds=bounds)
    ks = counts.keys()
    ks.sort()
    message('Result Counts in test set:')
    for k in ks:
      message(str((k, counts[k])))
  nExamples = len(trainExamples)
  message('Training with %d examples' % (nExamples))

  nVars = data.GetNVars()
  attrs = list(range(1, nVars + 1))
  nPossibleVals = data.GetNPossibleVals()
  for i in range(1, len(nPossibleVals)):
    if nPossibleVals[i - 1] == -1:
      attrs.remove(i)

  if details.pickleDataFileName != '':
    pickleDataFile = open(details.pickleDataFileName, 'wb+')
    pickle.dump(trainExamples, pickleDataFile)
    pickle.dump(testExamples, pickleDataFile)
    pickleDataFile.close()

  if details.bayesModel:
    composite = BayesComposite.BayesComposite()
  else:
    composite = Composite.Composite()

  composite._randomSeed = seed
  composite._splitFrac = details.splitFrac
  composite._shuffleActivities = details.shuffleActivities
  composite._randomizeActivities = details.randomActivities

  if hasattr(details, 'filterFrac'):
    composite._filterFrac = details.filterFrac
  if hasattr(details, 'filterVal'):
    composite._filterVal = details.filterVal

  composite.SetModelFilterData(details.modelFilterFrac, details.modelFilterVal)

  composite.SetActivityQuantBounds(details.activityBounds)
  nPossibleVals = data.GetNPossibleVals()
  if details.activityBounds:
    nPossibleVals[-1] = len(details.activityBounds) + 1

  if setDescNames:
    composite.SetInputOrder(data.GetVarNames())
    composite.SetDescriptorNames(details._descNames)
  else:
    composite.SetDescriptorNames(data.GetVarNames())
  composite.SetActivityQuantBounds(details.activityBounds)
  if details.nModels == 1:
    details.internalHoldoutFrac = 0.0
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

    composite.SetQuantBounds(details.qBounds)
    nPossibleVals = data.GetNPossibleVals()
    if details.activityBounds:
      nPossibleVals[-1] = len(details.activityBounds) + 1
    composite.Grow(
      trainExamples, attrs, nPossibleVals=[0] + nPossibleVals, buildDriver=driver, pruner=pruner,
      nTries=details.nModels, pruneIt=details.pruneIt, lessGreedy=details.lessGreedy,
      needsQuantization=0, treeBuilder=builder, nQuantBounds=details.qBounds,
      startAt=details.startAt, maxDepth=details.limitDepth, progressCallback=progressCallback,
      holdOutFrac=details.internalHoldoutFrac, replacementSelection=details.replacementSelection,
      recycleVars=details.recycleVars, randomDescriptors=details.randomDescriptors,
      silent=not _verbose)

  elif details.useSigTrees:
    from rdkit.ML.DecTree import CrossValidate
    from rdkit.ML.DecTree import BuildSigTree
    builder = BuildSigTree.SigTreeBuilder
    driver = CrossValidate.CrossValidationDriver
    nPossibleVals = data.GetNPossibleVals()
    if details.activityBounds:
      nPossibleVals[-1] = len(details.activityBounds) + 1
    if hasattr(details, 'sigTreeBiasList'):
      biasList = details.sigTreeBiasList
    else:
      biasList = None
    if hasattr(details, 'useCMIM'):
      useCMIM = details.useCMIM
    else:
      useCMIM = 0
    if hasattr(details, 'allowCollections'):
      allowCollections = details.allowCollections
    else:
      allowCollections = False
    composite.Grow(
      trainExamples, attrs, nPossibleVals=[0] + nPossibleVals, buildDriver=driver,
      nTries=details.nModels, needsQuantization=0, treeBuilder=builder, maxDepth=details.limitDepth,
      progressCallback=progressCallback, holdOutFrac=details.internalHoldoutFrac,
      replacementSelection=details.replacementSelection, recycleVars=details.recycleVars,
      randomDescriptors=details.randomDescriptors, biasList=biasList, useCMIM=useCMIM,
      allowCollection=allowCollections, silent=not _verbose)

  elif details.useKNN:
    from rdkit.ML.KNN import CrossValidate
    from rdkit.ML.KNN import DistFunctions

    driver = CrossValidate.CrossValidationDriver
    dfunc = ''
    if (details.knnDistFunc == "Euclidean"):
      dfunc = DistFunctions.EuclideanDist
    elif (details.knnDistFunc == "Tanimoto"):
      dfunc = DistFunctions.TanimotoDist
    else:
      assert 0, "Bad KNN distance metric value"

    composite.Grow(trainExamples, attrs, nPossibleVals=[0] + nPossibleVals, buildDriver=driver,
                   nTries=details.nModels, needsQuantization=0, numNeigh=details.knnNeighs,
                   holdOutFrac=details.internalHoldoutFrac, distFunc=dfunc)

  elif details.useNaiveBayes or details.useSigBayes:
    from rdkit.ML.NaiveBayes import CrossValidate
    driver = CrossValidate.CrossValidationDriver
    if not (hasattr(details, 'useSigBayes') and details.useSigBayes):
      composite.Grow(trainExamples, attrs, nPossibleVals=[0] + nPossibleVals, buildDriver=driver,
                     nTries=details.nModels, needsQuantization=0, nQuantBounds=details.qBounds,
                     holdOutFrac=details.internalHoldoutFrac,
                     replacementSelection=details.replacementSelection,
                     mEstimateVal=details.mEstimateVal, silent=not _verbose)
    else:
      if hasattr(details, 'useCMIM'):
        useCMIM = details.useCMIM
      else:
        useCMIM = 0

      composite.Grow(trainExamples, attrs, nPossibleVals=[0] + nPossibleVals, buildDriver=driver,
                     nTries=details.nModels, needsQuantization=0, nQuantBounds=details.qBounds,
                     mEstimateVal=details.mEstimateVal, useSigs=True, useCMIM=useCMIM,
                     holdOutFrac=details.internalHoldoutFrac,
                     replacementSelection=details.replacementSelection, silent=not _verbose)

    # #   elif details.useSVM:
    # #     from rdkit.ML.SVM import CrossValidate
    # #     driver = CrossValidate.CrossValidationDriver
    # #     composite.Grow(trainExamples, attrs, nPossibleVals=[0]+nPossibleVals,
    # #                    buildDriver=driver, nTries=details.nModels,
    # #                    needsQuantization=0,
    # #                    cost=details.svmCost,gamma=details.svmGamma,
    # #                    weights=details.svmWeights,degree=details.svmDegree,
    # #                    type=details.svmType,kernelType=details.svmKernel,
    # #                    coef0=details.svmCoeff,eps=details.svmEps,nu=details.svmNu,
    # #                    cache_size=details.svmCache,shrinking=details.svmShrink,
    # #                    dataType=details.svmDataType,
    # #                    holdOutFrac=details.internalHoldoutFrac,
    # #                    replacementSelection=details.replacementSelection,
    # #                    silent=not _verbose)

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
  message('# Overall Average Error: %%% 5.2f, Average Deviation: %%% 6.2f' %
          (100. * averageErr, 100. * avgDev))

  if details.bayesModel:
    composite.Train(trainExamples, verbose=0)

  # blow out the saved examples and then save the composite:
  composite.ClearModelExamples()
  if saveIt:
    composite.Pickle(details.outName)
  details.model = DbModule.binaryHolder(pickle.dumps(composite))

  badExamples = []
  if not details.detailedRes and (not hasattr(details, 'noScreen') or not details.noScreen):
    if details.splitRun:
      message('Testing all hold-out examples')
      wrong = testall(composite, testExamples, badExamples)
      message('%d examples (%% %5.2f) were misclassified' % (len(wrong), 100. * float(len(wrong)) /
                                                             float(len(testExamples))))
      _runDetails.holdout_error = float(len(wrong)) / len(testExamples)
    else:
      message('Testing all examples')
      wrong = testall(composite, namedExamples, badExamples)
      message('%d examples (%% %5.2f) were misclassified' % (len(wrong), 100. * float(len(wrong)) /
                                                             float(len(namedExamples))))
      _runDetails.overall_error = float(len(wrong)) / len(namedExamples)

  if details.detailedRes:
    message('\nEntire data set:')
    resTup = ScreenComposite.ShowVoteResults(
      range(data.GetNPts()), data, composite, nPossibleVals[-1], details.threshold)
    nGood, nBad, nSkip, avgGood, avgBad, avgSkip, voteTab = resTup
    nPts = len(namedExamples)
    nClass = nGood + nBad
    _runDetails.overall_error = float(nBad) / nClass
    _runDetails.overall_correct_conf = avgGood
    _runDetails.overall_incorrect_conf = avgBad
    _runDetails.overall_result_matrix = repr(voteTab)
    nRej = nClass - nPts
    if nRej > 0:
      _runDetails.overall_fraction_dropped = float(nRej) / nPts

    if details.splitRun:
      message('\nHold-out data:')
      resTup = ScreenComposite.ShowVoteResults(
        range(len(testExamples)), testExamples, composite, nPossibleVals[-1], details.threshold)
      nGood, nBad, nSkip, avgGood, avgBad, avgSkip, voteTab = resTup
      nPts = len(testExamples)
      nClass = nGood + nBad
      _runDetails.holdout_error = float(nBad) / nClass
      _runDetails.holdout_correct_conf = avgGood
      _runDetails.holdout_incorrect_conf = avgBad
      _runDetails.holdout_result_matrix = repr(voteTab)
      nRej = nClass - nPts
      if nRej > 0:
        _runDetails.holdout_fraction_dropped = float(nRej) / nPts

  if details.persistTblName and details.dbName:
    message('Updating results table %s:%s' % (details.dbName, details.persistTblName))
    details.Store(db=details.dbName, table=details.persistTblName)

  if details.badName != '':
    badFile = open(details.badName, 'w+')
    for i in range(len(badExamples)):
      ex = badExamples[i]
      vote = wrong[i]
      outStr = '%s\t%s\n' % (ex, vote)
      badFile.write(outStr)
    badFile.close()

  composite.ClearModelExamples()
  return composite


def RunIt(details, progressCallback=None, saveIt=1, setDescNames=0):
  """ does the actual work of building a composite model

    **Arguments**

      - details:  a _CompositeRun.CompositeRun_ object containing details
        (options, parameters, etc.) about the run

      - progressCallback: (optional) a function which is called with a single
        argument (the number of models built so far) after each model is built.

      - saveIt: (optional) if this is nonzero, the resulting model will be pickled
        and dumped to the filename specified in _details.outName_

      - setDescNames: (optional) if nonzero, the composite's _SetInputOrder()_ method
        will be called using the results of the data set's _GetVarNames()_ method;
        it is assumed that the details object has a _descNames attribute which
        is passed to the composites _SetDescriptorNames()_ method.  Otherwise
        (the default), _SetDescriptorNames()_ gets the results of _GetVarNames()_.

    **Returns**

      the composite model constructed


  """
  details.rundate = time.asctime()

  fName = details.tableName.strip()
  if details.outName == '':
    details.outName = fName + '.pkl'
  if not details.dbName:
    if details.qBounds != []:
      data = DataUtils.TextFileToData(fName)
    else:
      data = DataUtils.BuildQuantDataSet(fName)
  elif details.useSigTrees or details.useSigBayes:
    details.tableName = fName
    data = details.GetDataSet(pickleCol=0, pickleClass=DataStructs.ExplicitBitVect)
  elif details.qBounds != [] or not details.useTrees:
    details.tableName = fName
    data = details.GetDataSet()
  else:
    data = DataUtils.DBToQuantData(details.dbName,  # Function no longer defined
                                   fName,
                                   quantName=details.qTableName,
                                   user=details.dbUser,
                                   password=details.dbPassword)

  composite = RunOnData(details, data, progressCallback=progressCallback, saveIt=saveIt,
                        setDescNames=setDescNames)
  return composite


def ShowVersion(includeArgs=0):
  """ prints the version number

  """
  print('This is BuildComposite.py version %s' % (__VERSION_STRING))
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
  args, extra = getopt.getopt(
    sys.argv[1:],
    'P:o:n:p:b:sf:F:v:hlgd:rSTt:BQ:q:DVG:N:L:',
    ['nRuns=',
     'prune',
     'profile',
     'seed=',
     'noScreen',
     'modelFiltFrac=',
     'modelFiltVal=',
     'recycle',
     'randomDescriptors=',
     'doKnn',
     'knnK=',
     'knnTanimoto',
     'knnEuclid',
     'doSigTree',
     'allowCollections',
     'doNaiveBayes',
     'mEstimateVal=',
     'doSigBayes',

     # #                               'doSVM','svmKernel=','svmType=','svmGamma=',
     # #                               'svmCost=','svmWeights=','svmDegree=',
     # #                               'svmCoeff=','svmEps=','svmNu=','svmCache=',
     # #                               'svmShrink','svmDataType=',
     'replacementSelection', ])
  runDetails.profileIt = 0
  for arg, val in args:
    if arg == '-n':
      runDetails.nModels = int(val)
    elif arg == '-N':
      runDetails.note = val
    elif arg == '-o':
      runDetails.outName = val
    elif arg == '-Q':
      qBounds = eval(val)
      assert type(qBounds) in [type([]), type(
        ())], 'bad argument type for -Q, specify a list as a string'
      runDetails.activityBounds = qBounds
      runDetails.activityBoundsVals = val
    elif arg == '-p':
      runDetails.persistTblName = val
    elif arg == '-P':
      runDetails.pickleDataFileName = val
    elif arg == '-r':
      runDetails.randomActivities = 1
    elif arg == '-S':
      runDetails.shuffleActivities = 1
    elif arg == '-b':
      runDetails.badName = val
    elif arg == '-B':
      runDetails.bayesModels = 1
    elif arg == '-s':
      runDetails.splitRun = 1
    elif arg == '-f':
      runDetails.splitFrac = float(val)
    elif arg == '-F':
      runDetails.filterFrac = float(val)
    elif arg == '-v':
      runDetails.filterVal = float(val)
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
      assert type(qBounds) in [type([]), type(
        ())], 'bad argument type for -q, specify a list as a string'
      runDetails.qBoundCount = val
      runDetails.qBounds = qBounds
    elif arg == '-V':
      ShowVersion()
      sys.exit(0)
    elif arg == '--nRuns':
      runDetails.nRuns = int(val)
    elif arg == '--modelFiltFrac':
      runDetails.modelFilterFrac = float(val)
    elif arg == '--modelFiltVal':
      runDetails.modelFilterVal = float(val)
    elif arg == '--prune':
      runDetails.pruneIt = 1
    elif arg == '--profile':
      runDetails.profileIt = 1

    elif arg == '--recycle':
      runDetails.recycleVars = 1
    elif arg == '--randomDescriptors':
      runDetails.randomDescriptors = int(val)

    elif arg == '--doKnn':
      runDetails.useKNN = 1
      runDetails.useTrees = 0
      # #      runDetails.useSVM=0
      runDetails.useNaiveBayes = 0
    elif arg == '--knnK':
      runDetails.knnNeighs = int(val)
    elif arg == '--knnTanimoto':
      runDetails.knnDistFunc = "Tanimoto"
    elif arg == '--knnEuclid':
      runDetails.knnDistFunc = "Euclidean"

    elif arg == '--doSigTree':
      # #      runDetails.useSVM=0
      runDetails.useKNN = 0
      runDetails.useTrees = 0
      runDetails.useNaiveBayes = 0
      runDetails.useSigTrees = 1
    elif arg == '--allowCollections':
      runDetails.allowCollections = True

    elif arg == '--doNaiveBayes':
      runDetails.useNaiveBayes = 1
      # #      runDetails.useSVM=0
      runDetails.useKNN = 0
      runDetails.useTrees = 0
      runDetails.useSigBayes = 0
    elif arg == '--doSigBayes':
      runDetails.useSigBayes = 1
      runDetails.useNaiveBayes = 0
      # #      runDetails.useSVM=0
      runDetails.useKNN = 0
      runDetails.useTrees = 0
    elif arg == '--mEstimateVal':
      runDetails.mEstimateVal = float(val)

# #     elif arg == '--doSVM':
# #       runDetails.useSVM=1
# #       runDetails.useKNN=0
# #       runDetails.useTrees=0
# #       runDetails.useNaiveBayes=0
# #     elif arg == '--svmKernel':
# #       if val not in SVM.kernels.keys():
# #         message('kernel %s not in list of available kernels:\n%s\n'%(val,SVM.kernels.keys()))
# #         sys.exit(-1)
# #       else:
# #         runDetails.svmKernel=SVM.kernels[val]
# #     elif arg == '--svmType':
# #       if val not in SVM.machineTypes.keys():
# #         message('type %s not in list of available machines:\n%s\n'%(val,
# #                                                                     SVM.machineTypes.keys()))
# #         sys.exit(-1)
# #       else:
# #         runDetails.svmType=SVM.machineTypes[val]
# #     elif arg == '--svmGamma':
# #       runDetails.svmGamma = float(val)
# #     elif arg == '--svmCost':
# #       runDetails.svmCost = float(val)
# #     elif arg == '--svmWeights':
# #       # FIX: this is dangerous
# #       runDetails.svmWeights = eval(val)
# #     elif arg == '--svmDegree':
# #       runDetails.svmDegree = int(val)
# #     elif arg == '--svmCoeff':
# #       runDetails.svmCoeff = float(val)
# #     elif arg == '--svmEps':
# #       runDetails.svmEps = float(val)
# #     elif arg == '--svmNu':
# #       runDetails.svmNu = float(val)
# #     elif arg == '--svmCache':
# #       runDetails.svmCache = int(val)
# #     elif arg == '--svmShrink':
# #       runDetails.svmShrink = 0
# #     elif arg == '--svmDataType':
# #       runDetails.svmDataType=val

    elif arg == '--seed':
      # FIX: dangerous
      runDetails.randomSeed = eval(val)

    elif arg == '--noScreen':
      runDetails.noScreen = 1

    elif arg == '--replacementSelection':
      runDetails.replacementSelection = 1

    elif arg == '-h':
      Usage()

    else:
      Usage()
  runDetails.tableName = extra[0]

if __name__ == '__main__':
  if len(sys.argv) < 2:
    Usage()

  _runDetails.cmd = ' '.join(sys.argv)
  SetDefaults(_runDetails)
  ParseArgs(_runDetails)

  ShowVersion(includeArgs=1)

  if _runDetails.nRuns > 1:
    for i in range(_runDetails.nRuns):
      sys.stderr.write(
        '---------------------------------\n\tDoing %d of %d\n---------------------------------\n' %
        (i + 1, _runDetails.nRuns))
      RunIt(_runDetails)
  else:
    if _runDetails.profileIt:
      try:
        import hotshot
        import hotshot.stats
        prof = hotshot.Profile('prof.dat')
        prof.runcall(RunIt, _runDetails)
        stats = hotshot.stats.load('prof.dat')
        stats.strip_dirs()
        stats.sort_stats('time', 'calls')
        stats.print_stats(30)
      except ImportError:
        print('Profiling requires the hotshot module')
    else:
      RunIt(_runDetails)
