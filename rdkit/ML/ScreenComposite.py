# $Id$
#
#  Copyright (C) 2000-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" command line utility for screening composite models

**Usage**

  _ScreenComposite [optional args] modelfile(s) datafile_

Unless indicated otherwise (via command line arguments), _modelfile_ is
a file containing a pickled composite model and _filename_ is a QDAT file.

**Command Line Arguments**

  - -t *threshold value(s)*: use high-confidence predictions for the final
     analysis of the hold-out data.  The threshold value can be either a single
     float or a list/tuple of floats.  All thresholds should be between
     0.0 and 1.0

  - -D: do a detailed screen.

  - -d *database name*: instead of reading the data from a QDAT file,
     pull it from a database.  In this case, the _datafile_ argument
     provides the name of the database table containing the data set.

  - -N *note*: use all models from the database which have this note.
               The modelfile argument should contain the name of the table
               with the models.

  - -H: screen only the hold out set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -T: screen only the training set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -E: do a detailed Error analysis.  This shows each misclassified
     point and the number of times it was missed across all screened
     composites.  If the --enrich argument is also provided, only compounds
     that have true activity value equal to the enrichment value will be
     used.

  - --enrich *enrichVal*: target "active" value to be used in calculating
     enrichments.

  - -A: show All predictions.

  - -S: shuffle activity values before screening

  - -R: randomize activity values before screening

  - -F *filter frac*: filters the data before training to change the
     distribution of activity values in the training set.  *filter frac*
     is the fraction of the training set that should have the target value.
     **See note in BuildComposite help about data filtering**

  - -v *filter value*: filters the data before training to change the
     distribution of activity values in the training set. *filter value*
     is the target value to use in filtering.
     **See note in BuildComposite help about data filtering**

  - -V: be verbose when screening multiple models

  - -h: show this message and exit

  - --OOB: Do out an "out-of-bag" generalization error estimate.  This only
      makes sense when applied to the original data set.

  - --pickleCol *colId*: index of the column containing a pickled value
      (used primarily for cases where fingerprints are used as descriptors)

  *** Options for making Prediction (Hanneke) Plots ***

  - --predPlot=<fileName>: triggers the generation of a Hanneke plot and
      sets the name of the .txt file which will hold the output data.
      A Gnuplot control file, <fileName>.gnu, will also be generated.

  - --predActTable=<name> (optional):  name of the database table
      containing activity values.  If this is not provided, activities
      will be read from the same table containing the screening data

  - --predActCol=<name> (optional):  name of the activity column. If not
      provided, the name of the last column in the activity table will
      be used.

  - --predLogScale (optional):  If provided, the x axis of the
      prediction plot (the activity axis) will be plotted using a log
      scale

  - --predShow: launch a gnuplot instance and display the prediction
      plot (the plot will still be written to disk).

  *** The following options are likely obsolete ***

  - -P: read pickled data.  The datafile argument should contain
     a pickled data set. *relevant only to qdat files*

  - -q: data are not quantized (the composite should take care of
     quantization itself if it requires quantized data). *relevant only to
     qdat files*



"""


import os
import sys

import numpy

from rdkit import DataStructs
from rdkit.Dbase import DbModule
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import CompositeRun
from rdkit.ML.Data import DataUtils, SplitData
import pickle


try:
  from PIL import Image, ImageDraw
except ImportError:
  hasPil = 0
else:
  hasPil = 1

_details = CompositeRun.CompositeRun()

__VERSION_STRING = "3.3.0"


def message(msg, noRet=0):
  """ emits messages to _sys.stdout_
    override this in modules which import this one to redirect output

    **Arguments**

      - msg: the string to be displayed

  """
  if noRet:
    sys.stdout.write('%s ' % (msg))
  else:
    sys.stdout.write('%s\n' % (msg))


def error(msg):
  """ emits messages to _sys.stderr_
    override this in modules which import this one to redirect output

    **Arguments**

      - msg: the string to be displayed

  """
  sys.stderr.write('ERROR: %s\n' % (msg))


def CalcEnrichment(mat, tgt=1):
  if tgt < 0 or tgt >= mat.shape[0]:
    return 0
  nPts = float(sum(sum(mat)))
  nTgtPred = float(sum(mat[:, tgt]))
  if nTgtPred:
    pctCorrect = mat[tgt, tgt] / nTgtPred
    nTgtReal = float(sum(mat[tgt, :]))
    pctOverall = nTgtReal / nPts
  else:
    return 0.0
  return pctCorrect / pctOverall


def CollectResults(indices, dataSet, composite, callback=None, appendExamples=0, errorEstimate=0):
  """ screens a set of examples through a composite and returns the
  results
#DOC

  **Arguments**

    - examples: the examples to be screened (a sequence of sequences)
       it's assumed that the last element in each example is it's "value"

    - composite:  the composite model to be used

    - callback: (optional)  if provided, this should be a function
      taking a single argument that is called after each example is
      screened with the number of examples screened so far as the
      argument.

    - appendExamples: (optional)  this value is passed on to the
      composite's _ClassifyExample()_ method.

    - errorEstimate: (optional) calculate the "out of bag" error
      estimate for the composite using Breiman's definition.  This
      only makes sense when screening the original data set!
      [L. Breiman "Out-of-bag Estimation", UC Berkeley Dept of
      Statistics Technical Report (1996)]

  **Returns**

    a list of 3-tuples _nExamples_ long:

      1)  answer: the value from the example

      2)  pred: the composite model's prediction

      3)  conf: the confidence of the composite

  """
  # for i in range(len(composite)):
  #  print('  ',i,'TRAIN:',composite[i][0]._trainIndices)

  for j in range(len(composite)):
    tmp = composite.GetModel(j)
    if hasattr(tmp, '_trainIndices') and type(tmp._trainIndices) != dict:
      tis = {}
      if hasattr(tmp, '_trainIndices'):
        for v in tmp._trainIndices:
          tis[v] = 1
      tmp._trainIndices = tis

  nPts = len(indices)
  res = [None] * nPts
  for i in range(nPts):
    idx = indices[i]
    example = dataSet[idx]
    if errorEstimate:
      use = []
      for j in range(len(composite)):
        mdl = composite.GetModel(j)
        if not mdl._trainIndices.get(idx, 0):
          use.append(j)
    else:
      use = None
    # print('IDX:',idx,'use:',use  )
    pred, conf = composite.ClassifyExample(example, appendExample=appendExamples, onlyModels=use)
    if composite.GetActivityQuantBounds():
      answer = composite.QuantizeActivity(example)[-1]
    else:
      answer = example[-1]
    res[i] = answer, pred, conf
    if callback:
      callback(i)
  return res


def DetailedScreen(indices, data, composite, threshold=0, screenResults=None, goodVotes=None,
                   badVotes=None, noVotes=None, callback=None, appendExamples=0, errorEstimate=0):
  """ screens a set of examples cross a composite and breaks the
      predictions into *correct*,*incorrect* and *unclassified* sets.
#DOC
  **Arguments**

    - examples: the examples to be screened (a sequence of sequences)
       it's assumed that the last element in each example is its "value"

    - composite:  the composite model to be used

    - threshold: (optional) the threshold to be used to decide whether
      or not a given prediction should be kept

    - screenResults: (optional) the results of screening the results
      (a sequence of 3-tuples in the format returned by
      _CollectResults()_).  If this is provided, the examples will not
      be screened again.

    - goodVotes,badVotes,noVotes: (optional)  if provided these should
      be lists (or anything supporting an _append()_ method) which
      will be used to pass the screening results back.

    - callback: (optional)  if provided, this should be a function
      taking a single argument that is called after each example is
      screened with the number of examples screened so far as the
      argument.

    - appendExamples: (optional)  this value is passed on to the
      composite's _ClassifyExample()_ method.

    - errorEstimate: (optional) calculate the "out of bag" error
      estimate for the composite using Breiman's definition.  This
      only makes sense when screening the original data set!
      [L. Breiman "Out-of-bag Estimation", UC Berkeley Dept of
      Statistics Technical Report (1996)]

  **Notes**

    - since this function doesn't return anything, if one or more of
      the arguments _goodVotes_, _badVotes_, and _noVotes_ is not
      provided, there's not much reason to call it

  """
  if screenResults is None:
    screenResults = CollectResults(indices, data, composite, callback=callback,
                                   appendExamples=appendExamples, errorEstimate=errorEstimate)
  if goodVotes is None:
    goodVotes = []
  if badVotes is None:
    badVotes = []
  if noVotes is None:
    noVotes = []
  for i in range(len(screenResults)):
    answer, pred, conf = screenResults[i]
    if conf > threshold:
      if pred != answer:
        badVotes.append((answer, pred, conf, i))
      else:
        goodVotes.append((answer, pred, conf, i))
    else:
      noVotes.append((answer, pred, conf, i))


def ShowVoteResults(indices, data, composite, nResultCodes, threshold, verbose=1,
                    screenResults=None, callback=None, appendExamples=0, goodVotes=None,
                    badVotes=None, noVotes=None, errorEstimate=0):
  """ screens the results and shows a detailed workup

  The work of doing the screening and processing the results is
  handled by _DetailedScreen()_
#DOC

  **Arguments**

    - examples: the examples to be screened (a sequence of sequences)
       it's assumed that the last element in each example is its "value"

    - composite:  the composite model to be used

    - nResultCodes: the number of possible results the composite can
      return

    - threshold: the threshold to be used to decide whether or not a
      given prediction should be kept

    - screenResults: (optional) the results of screening the results
      (a sequence of 3-tuples in the format returned by
      _CollectResults()_).  If this is provided, the examples will not
      be screened again.

    - callback: (optional)  if provided, this should be a function
      taking a single argument that is called after each example is
      screened with the number of examples screened so far as the
      argument.

    - appendExamples: (optional)  this value is passed on to the
      composite's _ClassifyExample()_ method.

    - goodVotes,badVotes,noVotes: (optional)  if provided these should
      be lists (or anything supporting an _append()_ method) which
      will be used to pass the screening results back.

    - errorEstimate: (optional) calculate the "out of bag" error
      estimate for the composite using Breiman's definition.  This
      only makes sense when screening the original data set!
      [L. Breiman "Out-of-bag Estimation", UC Berkeley Dept of
      Statistics Technical Report (1996)]

  **Returns**

    a 7-tuple:

      1) the number of good (correct) predictions

      2) the number of bad (incorrect) predictions

      3) the number of predictions skipped due to the _threshold_

      4) the average confidence in the good predictions

      5) the average confidence in the bad predictions

      6) the average confidence in the skipped predictions

      7) the results table

  """
  nExamples = len(indices)
  if goodVotes is None:
    goodVotes = []
  if badVotes is None:
    badVotes = []
  if noVotes is None:
    noVotes = []
  DetailedScreen(indices, data, composite, threshold, screenResults=screenResults,
                 goodVotes=goodVotes, badVotes=badVotes, noVotes=noVotes, callback=callback,
                 appendExamples=appendExamples, errorEstimate=errorEstimate)
  nBad = len(badVotes)
  nGood = len(goodVotes)
  nClassified = nGood + nBad
  if verbose:
    print('\n\t*** Vote Results ***')
    print('misclassified: %d/%d (%%%4.2f)\t%d/%d (%%%4.2f)' %
          (nBad, nExamples, 100. * float(nBad) / nExamples, nBad, nClassified,
           100. * float(nBad) / nClassified))
  nSkip = len(noVotes)
  if nSkip > 0:
    if verbose:
      print('skipped: %d/%d (%%% 4.2f)' % (nSkip, nExamples, 100. * float(nSkip) / nExamples))
    noConf = numpy.array([x[2] for x in noVotes])
    avgSkip = sum(noConf) / float(nSkip)
  else:
    avgSkip = 0.

  if nBad > 0:
    badConf = numpy.array([x[2] for x in badVotes])
    avgBad = sum(badConf) / float(nBad)
  else:
    avgBad = 0.

  if nGood > 0:
    goodRes = [x[1] for x in goodVotes]
    goodConf = numpy.array([x[2] for x in goodVotes])
    avgGood = sum(goodConf) / float(nGood)
  else:
    goodRes = []
    goodConf = []
    avgGood = 0.

  if verbose:
    print()
    print('average correct confidence:   % 6.4f' % avgGood)
    print('average incorrect confidence: % 6.4f' % avgBad)

  voteTab = numpy.zeros((nResultCodes, nResultCodes), numpy.int)
  for res in goodRes:
    voteTab[res, res] += 1
  for ans, res, conf, idx in badVotes:
    voteTab[ans, res] += 1

  if verbose:
    print()
    print('\tResults Table:')
    vTab = voteTab.transpose()
    colCounts = numpy.sum(vTab, 0)
    rowCounts = numpy.sum(vTab, 1)
    message('')
    for i in range(nResultCodes):
      if rowCounts[i] == 0:
        rowCounts[i] = 1
      row = vTab[i]
      message('    ', noRet=1)
      for j in range(nResultCodes):
        entry = row[j]
        message(' % 6d' % entry, noRet=1)
      message('     | % 4.2f' % (100. * vTab[i, i] / rowCounts[i]))
    message('    ', noRet=1)
    for i in range(nResultCodes):
      message('-------', noRet=1)
    message('')
    message('    ', noRet=1)
    for i in range(nResultCodes):
      if colCounts[i] == 0:
        colCounts[i] = 1
      message(' % 6.2f' % (100. * vTab[i, i] / colCounts[i]), noRet=1)
    message('')

  return nGood, nBad, nSkip, avgGood, avgBad, avgSkip, voteTab


def ScreenIt(composite, indices, data, partialVote=0, voteTol=0.0, verbose=1, screenResults=None,
             goodVotes=None, badVotes=None, noVotes=None):
  """ screens a set of data using a composite model and prints out
             statistics about the screen.
#DOC
    The work of doing the screening and processing the results is
    handled by _DetailedScreen()_

  **Arguments**

    - composite:  the composite model to be used

    - data: the examples to be screened (a sequence of sequences)
       it's assumed that the last element in each example is its "value"

    - partialVote: (optional) toggles use of the threshold value in
      the screnning.

    - voteTol: (optional) the threshold to be used to decide whether or not a
      given prediction should be kept

    - verbose: (optional) sets degree of verbosity of the screening

    - screenResults: (optional) the results of screening the results
      (a sequence of 3-tuples in the format returned by
      _CollectResults()_).  If this is provided, the examples will not
      be screened again.

    - goodVotes,badVotes,noVotes: (optional)  if provided these should
      be lists (or anything supporting an _append()_ method) which
      will be used to pass the screening results back.


  **Returns**

    a 7-tuple:

      1) the number of good (correct) predictions

      2) the number of bad (incorrect) predictions

      3) the number of predictions skipped due to the _threshold_

      4) the average confidence in the good predictions

      5) the average confidence in the bad predictions

      6) the average confidence in the skipped predictions

      7) None

  """
  if goodVotes is None:
    goodVotes = []
  if badVotes is None:
    badVotes = []
  if noVotes is None:
    noVotes = []

  if not partialVote:
    voteTol = 0.0

  DetailedScreen(indices, data, composite, voteTol, screenResults=screenResults,
                 goodVotes=goodVotes, badVotes=badVotes, noVotes=noVotes)

  nGood = len(goodVotes)
  goodAccum = 0.
  for res, pred, conf, idx in goodVotes:
    goodAccum += conf

  misCount = len(badVotes)
  badAccum = 0.
  for res, pred, conf, idx in badVotes:
    badAccum += conf

  nSkipped = len(noVotes)
  goodSkipped = 0
  badSkipped = 0
  skipAccum = 0.
  for ans, pred, conf, idx in noVotes:
    skipAccum += conf
    if ans != pred:
      badSkipped += 1
    else:
      goodSkipped += 1

  nData = nGood + misCount + nSkipped
  if verbose:
    print('Total N Points:', nData)
  if partialVote:
    nCounted = nData - nSkipped
    if verbose:
      print('Misclassifications: %d (%%%4.2f)' % (misCount, 100. * float(misCount) / nCounted))
      print('N Skipped: %d (%%%4.2f)' % (nSkipped, 100. * float(nSkipped) / nData))
      print('\tGood Votes Skipped: %d (%%%4.2f)' %
            (goodSkipped, 100. * float(goodSkipped) / nSkipped))
      print('\tBad Votes Skipped: %d (%%%4.2f)' % (badSkipped, 100. * float(badSkipped) / nSkipped))
  else:
    if verbose:
      print('Misclassifications: %d (%%%4.2f)' % (misCount, 100. * float(misCount) / nData))
      print('Average Correct Vote Confidence:   % 6.4f' % (goodAccum / (nData - misCount)))
      print('Average InCorrect Vote Confidence: % 6.4f' % (badAccum / misCount))

  avgGood = 0
  avgBad = 0
  avgSkip = 0
  if nGood:
    avgGood = goodAccum / nGood
  if misCount:
    avgBad = badAccum / misCount
  if nSkipped:
    avgSkip = skipAccum / nSkipped
  return nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, None


def _processVoteList(votes, data):
  """ *Internal Use Only*

  converts a list of 4 tuples: (answer,prediction,confidence,idx) into
  an alternate list: (answer,prediction,confidence,data point)

   **Arguments**

     - votes: a list of 4 tuples: (answer, prediction, confidence,
       index)

     - data: a _DataUtils.MLData.MLDataSet_


   **Note**: alterations are done in place in the _votes_ list

  """
  for i in range(len(votes)):
    ans, pred, conf, idx = votes[i]
    votes[i] = (ans, pred, conf, data[idx])


def PrepareDataFromDetails(model, details, data, verbose=0):
  if (hasattr(details, 'doHoldout') and details.doHoldout) or \
     (hasattr(details, 'doTraining') and details.doTraining):
    try:
      splitF = model._splitFrac
    except AttributeError:
      pass
    else:
      if verbose:
        message('s', noRet=1)

      if hasattr(details, 'errorEstimate') and details.errorEstimate and \
         hasattr(details, 'doHoldout') and details.doHoldout:
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
        message('******  WARNING: OOB screening should not be combined with doHoldout option.')
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
      trainIdx, testIdx = SplitData.SplitIndices(data.GetNPts(), splitF, silent=1)

    if hasattr(details, 'filterFrac') and details.filterFrac != 0.0:
      if verbose:
        message('f', noRet=1)
      trainFilt, temp = DataUtils.FilterData(data, details.filterVal, details.filterFrac, -1,
                                             indicesToUse=trainIdx, indicesOnly=1)
      testIdx += temp
      trainIdx = trainFilt
  elif hasattr(details, 'errorEstimate') and details.errorEstimate:
    # the OOB screening works by checking to see if a given index
    # is in the
    if hasattr(details, 'filterFrac') and details.filterFrac != 0.0:
      if verbose:
        message('f', noRet=1)
      testIdx, trainIdx = DataUtils.FilterData(data, details.filterVal, details.filterFrac, -1,
                                               indicesToUse=range(data.GetNPts()), indicesOnly=1)
      testIdx.extend(trainIdx)
    else:
      testIdx = list(range(data.GetNPts()))
    trainIdx = []
  else:
    testIdx = list(range(data.GetNPts()))
    trainIdx = []
  if hasattr(details, 'doTraining') and details.doTraining:
    testIdx, trainIdx = trainIdx, testIdx
  return trainIdx, testIdx


def ScreenFromDetails(models, details, callback=None, setup=None, appendExamples=0, goodVotes=None,
                      badVotes=None, noVotes=None, data=None, enrichments=None):
  """  Screens a set of data using a a _CompositeRun.CompositeRun_
       instance to provide parameters

# DOC

    The actual data to be used are extracted from the database and
    table specified in _details_

    Aside from dataset construction,  _ShowVoteResults()_ does most of
    the heavy lifting here.

  **Arguments**

      - model: a composite model

      - details:  a _CompositeRun.CompositeRun_ object containing details
        (options, parameters, etc.) about the run

      - callback: (optional)  if provided, this should be a function
        taking a single argument that is called after each example is
        screened with the number of examples screened so far as the
        argument.

      - setup: (optional) a function taking a single argument which is
        called at the start of screening with the number of points to
        be screened as the argument.

      - appendExamples: (optional)  this value is passed on to the
        composite's _ClassifyExample()_ method.

      - goodVotes,badVotes,noVotes: (optional)  if provided these should
        be lists (or anything supporting an _append()_ method) which
        will be used to pass the screening results back.


  **Returns**

    a 7-tuple:

      1) the number of good (correct) predictions

      2) the number of bad (incorrect) predictions

      3) the number of predictions skipped due to the _threshold_

      4) the average confidence in the good predictions

      5) the average confidence in the bad predictions

      6) the average confidence in the skipped predictions

      7) the results table

  """
  if data is None:
    if hasattr(details, 'pickleCol'):
      data = details.GetDataSet(pickleCol=details.pickleCol,
                                pickleClass=DataStructs.ExplicitBitVect)
    else:
      data = details.GetDataSet()
  if details.threshold > 0.0:
    details.partialVote = 1
  else:
    details.partialVote = 0

  if type(models) not in [list, tuple]:
    models = (models, )

  nModels = len(models)

  if setup is not None:
    setup(nModels * data.GetNPts())

  nGood = numpy.zeros(nModels, numpy.float)
  nBad = numpy.zeros(nModels, numpy.float)
  nSkip = numpy.zeros(nModels, numpy.float)
  confGood = numpy.zeros(nModels, numpy.float)
  confBad = numpy.zeros(nModels, numpy.float)
  confSkip = numpy.zeros(nModels, numpy.float)
  voteTab = None
  if goodVotes is None:
    goodVotes = []
  if badVotes is None:
    badVotes = []
  if noVotes is None:
    noVotes = []
  if enrichments is None:
    enrichments = [0.0] * nModels
  badVoteDict = {}
  noVoteDict = {}

  for i in range(nModels):
    if nModels > 1:
      goodVotes = []
      badVotes = []
      noVotes = []
    model = models[i]

    try:
      seed = model._randomSeed
    except AttributeError:
      pass
    else:
      DataUtils.InitRandomNumbers(seed)

    if (hasattr(details, 'shuffleActivities') and details.shuffleActivities) or \
       (hasattr(details, 'randomActivities') and details.randomActivities):
      if hasattr(details, 'shuffleActivities') and details.shuffleActivities:
        shuffle = True
      else:
        shuffle = False
      randomize = True
      DataUtils.RandomizeActivities(data, shuffle=shuffle, runDetails=details)
    else:
      randomize = False
      shuffle = False

    if hasattr(model, '_shuffleActivities') and \
       model._shuffleActivities and \
       not shuffle:
      message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
      message('******  WARNING: Shuffled model being screened with unshuffled data.')
      message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
    if hasattr(model, '_randomizeActivities') and \
       model._randomizeActivities and \
       not randomize:
      message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
      message('******  WARNING: Random model being screened with non-random data.')
      message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')

    trainIdx, testIdx = PrepareDataFromDetails(model, details, data)

    nPossible = model.GetQuantBounds()[1]
    if callback:
      cb = lambda x, y=callback, z=i * data.GetNPts(): y(x + z)
    else:
      cb = None
    if not hasattr(details, 'errorEstimate') or not details.errorEstimate:
      errorEstimate = 0
    else:
      errorEstimate = 1
    g, b, s, aG, aB, aS, vT = ShowVoteResults(
      testIdx, data, model, nPossible[-1], details.threshold, verbose=0, callback=cb,
      appendExamples=appendExamples, goodVotes=goodVotes, badVotes=badVotes, noVotes=noVotes,
      errorEstimate=errorEstimate)
    if voteTab is None:
      voteTab = numpy.zeros(vT.shape, numpy.float)
    if hasattr(details, 'errorAnalysis') and details.errorAnalysis:
      for a, p, c, idx in badVotes:
        label = testIdx[idx]
        if hasattr(details, 'enrichTgt') and details.enrichTgt >= 0:
          if a == details.enrichTgt:
            badVoteDict[label] = badVoteDict.get(label, 0) + 1
        else:
          badVoteDict[label] = badVoteDict.get(label, 0) + 1
      for a, p, c, idx in noVotes:
        label = testIdx[idx]
        if hasattr(details, 'enrichTgt') and details.enrichTgt >= 0:
          if a == details.enrichTgt:
            noVoteDict[label] = noVoteDict.get(label, 0) + 1
        else:
          noVoteDict[label] = noVoteDict.get(label, 0) + 1

    voteTab += vT
    nGood[i] = g
    nBad[i] = b
    nSkip[i] = s
    confGood[i] = aG
    confBad[i] = aB
    confSkip[i] = aS

    if hasattr(details, 'enrichTgt') and details.enrichTgt >= 0:
      enrichments[i] = CalcEnrichment(vT, tgt=details.enrichTgt)

  if nModels == 1:
    return g, b, s, aG, aB, aS, vT
  else:
    voteTab /= nModels

    avgNBad = sum(nBad) / nModels
    devNBad = numpy.sqrt(sum((nBad - avgNBad)**2) / (nModels - 1))

    # bestIdx = numpy.argsort(nBad)[0]

    avgNGood = sum(nGood) / nModels
    devNGood = numpy.sqrt(sum((nGood - avgNGood)**2) / (nModels - 1))

    avgNSkip = sum(nSkip) / nModels
    devNSkip = numpy.sqrt(sum((nSkip - avgNSkip)**2) / (nModels - 1))

    avgConfBad = sum(confBad) / nModels
    devConfBad = numpy.sqrt(sum((confBad - avgConfBad)**2) / (nModels - 1))

    avgConfGood = sum(confGood) / nModels
    devConfGood = numpy.sqrt(sum((confGood - avgConfGood)**2) / (nModels - 1))

    avgConfSkip = sum(confSkip) / nModels
    devConfSkip = numpy.sqrt(sum((confSkip - avgConfSkip)**2) / (nModels - 1))
    return ((avgNGood, devNGood), (avgNBad, devNBad), (avgNSkip, devNSkip),
            (avgConfGood, devConfGood), (avgConfBad, devConfBad), (avgConfSkip, devConfSkip),
            voteTab)


def GetScreenImage(nGood, nBad, nRej, size=None):
  if not hasPil:
    return None
  try:
    nTot = float(nGood) + float(nBad) + float(nRej)
  except TypeError:
    nGood = nGood[0]
    nBad = nBad[0]
    nRej = nRej[0]
    nTot = float(nGood) + float(nBad) + float(nRej)

  if not nTot:
    return None
  goodColor = (100, 100, 255)
  badColor = (255, 100, 100)
  rejColor = (255, 255, 100)

  pctGood = float(nGood) / nTot
  pctBad = float(nBad) / nTot
  pctRej = float(nRej) / nTot

  if size is None:
    size = (100, 100)
  img = Image.new('RGB', size, (255, 255, 255))
  draw = ImageDraw.Draw(img)
  box = (0, 0, size[0] - 1, size[1] - 1)

  startP = -90
  endP = int(startP + pctGood * 360)
  draw.pieslice(box, startP, endP, fill=goodColor)
  startP = endP
  endP = int(startP + pctBad * 360)
  draw.pieslice(box, startP, endP, fill=badColor)
  startP = endP
  endP = int(startP + pctRej * 360)
  draw.pieslice(box, startP, endP, fill=rejColor)

  return img


def ScreenToHtml(nGood, nBad, nRej, avgGood, avgBad, avgSkip, voteTable, imgDir='.', fullPage=1,
                 skipImg=0, includeDefs=1):
  """ returns the text of a web page showing the screening details
#DOC
    **Arguments**

     - nGood: number of correct predictions

     - nBad:  number of incorrect predictions

     - nRej:  number of rejected predictions

     - avgGood: average correct confidence

     - avgBad:  average incorrect confidence

     - avgSkip: average rejected confidence

     - voteTable: vote table

     - imgDir: (optional) the directory to be used to hold the vote
       image (if constructed)

   **Returns**

     a string containing HTML

  """
  if type(nGood) == tuple:
    multModels = 1
  else:
    multModels = 0

  if fullPage:
    outTxt = ["""<html><body>"""]
    outTxt.append('<center><h2>VOTE DETAILS</h2></center>')
  else:
    outTxt = []

  outTxt.append('<font>')

  # Get the image
  if not skipImg:
    img = GetScreenImage(nGood, nBad, nRej)
    if img:
      if imgDir:
        imgFileName = '/'.join((imgDir, 'votes.png'))
      else:
        imgFileName = 'votes.png'
      img.save(imgFileName)
      outTxt.append('<center><img src="%s"></center>' % (imgFileName))

  nPoss = len(voteTable)
  pureCounts = numpy.sum(voteTable, 1)
  accCounts = numpy.sum(voteTable, 0)
  pureVect = numpy.zeros(nPoss, numpy.float)
  accVect = numpy.zeros(nPoss, numpy.float)
  for i in range(nPoss):
    if pureCounts[i]:
      pureVect[i] = float(voteTable[i, i]) / pureCounts[i]
    if accCounts[i]:
      accVect[i] = float(voteTable[i, i]) / accCounts[i]

  outTxt.append('<center><table border=1>')
  outTxt.append('<tr><td></td>')
  for i in range(nPoss):
    outTxt.append('<th>%d</th>' % i)
  outTxt.append('<th>% Accurate</th>')
  outTxt.append('</tr>')
  # outTxt.append('<th rowspan=%d>Predicted</th></tr>'%(nPoss+1))
  for i in range(nPoss):
    outTxt.append('<tr><th>%d</th>' % (i))
    for j in range(nPoss):
      if i == j:
        if not multModels:
          outTxt.append('<td bgcolor="#A0A0FF">%d</td>' % (voteTable[j, i]))
        else:
          outTxt.append('<td bgcolor="#A0A0FF">%.2f</td>' % (voteTable[j, i]))
      else:
        if not multModels:
          outTxt.append('<td>%d</td>' % (voteTable[j, i]))
        else:
          outTxt.append('<td>%.2f</td>' % (voteTable[j, i]))
    outTxt.append('<td>%4.2f</td</tr>' % (100.0 * accVect[i]))
    if i == 0:
      outTxt.append('<th rowspan=%d>Predicted</th></tr>' % (nPoss))
    else:
      outTxt.append('</tr>')
  outTxt.append('<tr><th>% Pure</th>')
  for i in range(nPoss):
    outTxt.append('<td>%4.2f</td>' % (100.0 * pureVect[i]))
  outTxt.append('</tr>')
  outTxt.append('<tr><td></td><th colspan=%d>Original</th>' % (nPoss))
  outTxt.append('</table></center>')

  if not multModels:
    nTotal = nBad + nGood + nRej
    nClass = nBad + nGood
    if nClass:
      pctErr = 100. * float(nBad) / nClass
    else:
      pctErr = 0.0

    outTxt.append('<p>%d of %d examples were misclassified (%%%4.2f)' %
                  (nBad, nGood + nBad, pctErr))
    if nRej > 0:
      pctErr = 100. * float(nBad) / (nGood + nBad + nRej)
      outTxt.append('<p>                %d of %d overall: (%%%4.2f)' % (nBad, nTotal, pctErr))
      pctRej = 100. * float(nRej) / nTotal
      outTxt.append('<p>%d of %d examples were rejected (%%%4.2f)' % (nRej, nTotal, pctRej))
    if nGood != 0:
      outTxt.append('<p>The correctly classified examples had an average confidence of %6.4f' %
                    avgGood)

    if nBad != 0:
      outTxt.append('<p>The incorrectly classified examples had an average confidence of %6.4f' %
                    avgBad)
    if nRej != 0:
      outTxt.append('<p>The rejected examples had an average confidence of %6.4f' % avgSkip)
  else:
    nTotal = nBad[0] + nGood[0] + nRej[0]
    nClass = nBad[0] + nGood[0]
    devClass = nBad[1] + nGood[1]
    if nClass:
      pctErr = 100. * float(nBad[0]) / nClass
      devPctErr = 100. * float(nBad[1]) / nClass
    else:
      pctErr = 0.0
      devPctErr = 0.0

    outTxt.append('<p>%.2f(%.2f) of %.2f(%.2f) examples were misclassified (%%%4.2f(%4.2f))' %
                  (nBad[0], nBad[1], nClass, devClass, pctErr, devPctErr))
    if nRej > 0:
      pctErr = 100. * float(nBad[0]) / nTotal
      devPctErr = 100. * float(nBad[1]) / nTotal
      outTxt.append('<p>                %.2f(%.2f) of %d overall: (%%%4.2f(%4.2f))' %
                    (nBad[0], nBad[1], nTotal, pctErr, devPctErr))
      pctRej = 100. * float(nRej[0]) / nTotal
      devPctRej = 100. * float(nRej[1]) / nTotal
      outTxt.append('<p>%.2f(%.2f) of %d examples were rejected (%%%4.2f(%4.2f))' %
                    (nRej[0], nRej[1], nTotal, pctRej, devPctRej))
    if nGood != 0:
      outTxt.append(
        '<p>The correctly classified examples had an average confidence of %6.4f(%.4f)' % avgGood)

    if nBad != 0:
      outTxt.append(
        '<p>The incorrectly classified examples had an average confidence of %6.4f(%.4f)' % avgBad)
    if nRej != 0:
      outTxt.append('<p>The rejected examples had an average confidence of %6.4f(%.4f)' % avgSkip)

  outTxt.append('</font>')
  if includeDefs:
    txt = """
    <p><b>Definitions:</b>
    <ul>
    <li> <i>% Pure:</i>  The percentage of, for example, known positives predicted to be positive.
    <li> <i>% Accurate:</i>  The percentage of, for example, predicted positives that actually
      are positive.
    </ul>
    """
    outTxt.append(txt)

  if fullPage:
    outTxt.append("""</body></html>""")
  return '\n'.join(outTxt)


def MakePredPlot(details, indices, data, goodVotes, badVotes, nRes, idCol=0, verbose=0):
  """

  **Arguments**

    - details:  a CompositeRun.RunDetails object

    - indices: a sequence of integer indices into _data_

    - data: the data set in question.  We assume that the ids for
      the data points are in the _idCol_ column

    - goodVotes/badVotes: predictions where the model was correct/incorrect.
      These are sequences of 4-tuples:
        (answer,prediction,confidence,index into _indices_)

  """
  if not hasattr(details, 'predPlot') or not details.predPlot:
    return

  if verbose:
    message('\n-> Constructing Prediction (Hanneke) Plot')
  outF = open(details.predPlot, 'w+')
  gnuF = open('%s.gnu' % details.predPlot, 'w+')
  # first get the ids of the data points we screened:
  ptIds = [data[x][idCol] for x in indices]

  # get a connection to the database we'll use to grab the continuous
  #  activity values:
  origConn = DbConnect(details.dbName, details.tableName, user=details.dbUser,
                       password=details.dbPassword)
  colNames = origConn.GetColumnNames()
  idName = colNames[idCol]
  if not hasattr(details, 'predActTable') or \
     not details.predActTable or \
     details.predActTable == details.tableName:
    actConn = origConn
  else:
    actConn = DbConnect(details.dbName, details.predActTable, user=details.dbUser,
                        password=details.dbPassword)
  if verbose:
    message('\t-> Pulling Activity Data')

  if type(ptIds[0]) not in [type(''), type(u'')]:
    ptIds = [str(x) for x in ptIds]
  whereL = [DbModule.placeHolder] * len(ptIds)
  if hasattr(details, 'predActCol') and details.predActCol:
    actColName = details.predActCol
  else:
    actColName = actConn.GetColumnNames()[-1]

  whereTxt = "%s in (%s)" % (idName, ','.join(whereL))
  rawD = actConn.GetData(fields='%s,%s' % (idName, actColName), where=whereTxt, extras=ptIds)
  # order the data returned:
  if verbose:
    message('\t-> Creating Plot')
  acts = [None] * len(ptIds)
  for entry in rawD:
    ID, act = entry
    idx = ptIds.index(ID)
    acts[idx] = act
  outF.write('#ID Pred Conf %s\n' % (actColName))
  for ans, pred, conf, idx in goodVotes:
    act = acts[idx]
    if act != 'None':
      act = float(act)
    else:
      act = 0
    outF.write('%s %d %.4f %f\n' % (ptIds[idx], pred, conf, act))
  for ans, pred, conf, idx in badVotes:
    act = acts[idx]
    if act != 'None':
      act = float(act)
    else:
      act = 0
    outF.write('%s %d %.4f %f\n' % (ptIds[idx], pred, conf, act))
  outF.close()
  if not hasattr(details, 'predLogScale') or not details.predLogScale:
    actLabel = actColName
  else:
    actLabel = 'log(%s)' % (actColName)
  actLabel = actLabel.replace('_', ' ')
  gnuHdr = """# Generated by ScreenComposite.py version: %s
  set size square 0.7
  set yrange [:1]
  set data styl points
  set ylab 'confidence'
  set xlab '%s'
  set grid
  set nokey
  set term postscript enh color solid "Helvetica" 16
  set term X
  """ % (__VERSION_STRING, actLabel)
  gnuF.write(gnuHdr)
  plots = []
  for i in range(nRes):
    if not hasattr(details, 'predLogScale') or not details.predLogScale:
      plots.append("'%s' us 4:($2==%d?$3:0/0)" % (details.predPlot, i))
    else:
      plots.append("'%s' us (log10($4)):($2==%d?$3:0/0)" % (details.predPlot, i))
  gnuF.write("plot %s\n" % (','.join(plots)))
  gnuTail = """
  # EOF
  """
  gnuF.write(gnuTail)
  gnuF.close()
  if hasattr(details, 'predShow') and details.predShow:
    try:
      try:
        from Gnuplot import Gnuplot
      except ImportError:
        raise ImportError('Functionality requires the Gnuplot module')
      p = Gnuplot()
      p('cd "%s"' % (os.getcwd()))
      p('load "%s.gnu"' % (details.predPlot))
      input('press return to continue...\n')
    except Exception:
      import traceback
      traceback.print_exc()


def Go(details):
  pass


def SetDefaults(details=None):
  global _details
  if details is None:
    details = _details
  CompositeRun.SetDefaults(details)
  details.screenVoteTol = [0.]
  details.detailedScreen = 0
  details.doHoldout = 0
  details.doTraining = 0
  details.errorAnalysis = 0
  details.verbose = 0
  details.partialVote = 0
  return details


def Usage():
  """ prints a list of arguments for when this is used from the
  command line and then exits

  """
  print(__doc__)
  sys.exit(-1)


def ShowVersion(includeArgs=0):
  """ prints the version number of the program

  """
  print('This is ScreenComposite.py version %s' % (__VERSION_STRING))
  if includeArgs:
    print('command line was:')
    print(' '.join(sys.argv))


def ParseArgs(details):
  import getopt
  try:
    args, extras = getopt.getopt(sys.argv[1:], 'EDd:t:VN:HThSRF:v:AX', ['predPlot=',
                                                                        'predActCol=',
                                                                        'predActTable=',
                                                                        'predLogScale',
                                                                        'predShow',
                                                                        'OOB',
                                                                        'pickleCol=',
                                                                        'enrich=', ])
  except Exception:
    import traceback
    traceback.print_exc()
    Usage()

  details.predPlot = ''
  details.predActCol = ''
  details.predActTable = ''
  details.predLogScale = ''
  details.predShow = 0
  details.errorEstimate = 0
  details.pickleCol = -1
  details.enrichTgt = -1
  for arg, val in args:
    if arg == '-d':
      details.dbName = val
    elif arg == '-D':
      details.detailedScreen = 1
    elif arg == '-t':
      details.partialVote = 1
      voteTol = eval(val)
      if type(voteTol) not in [type([]), type((1, 1))]:
        voteTol = [voteTol]
      for tol in voteTol:
        if tol > 1 or tol < 0:
          error('Voting threshold must be between 0 and 1')
          sys.exit(-2)
      details.screenVoteTol = voteTol
    elif arg == '-N':
      details.note = val
    elif arg == '-H':
      details.doTraining = 0
      details.doHoldout = 1
    elif arg == '-T':
      details.doHoldout = 0
      details.doTraining = 1
    elif arg == '-E':
      details.errorAnalysis = 1
      details.detailedScreen = 1
    elif arg == '-A':
      details.showAll = 1
      details.detailedScreen = 1
    elif arg == '-S':
      details.shuffleActivities = 1
    elif arg == '-R':
      details.randomActivities = 1
    elif arg == '-h':
      Usage()
    elif arg == '-F':
      details.filterFrac = float(val)
    elif arg == '-v':
      details.filterVal = float(val)
    elif arg == '-V':
      verbose = 1
    elif arg == '--predPlot':
      details.detailedScreen = 1
      details.predPlot = val
    elif arg == '--predActCol':
      details.predActCol = val
    elif arg == '--predActTable':
      details.predActTable = val
    elif arg == '--predLogScale':
      details.predLogScale = 1
    elif arg == '--predShow':
      details.predShow = 1
    elif arg == '--predShow':
      details.predShow = 1
    elif arg == '--OOB':
      details.errorEstimate = 1
    elif arg == '--pickleCol':
      details.pickleCol = int(val) - 1
    elif arg == '--enrich':
      details.enrichTgt = int(val)
    else:
      Usage()

  if len(extras) < 1:
    Usage()
  return extras


if __name__ == '__main__':
  details = SetDefaults()
  extras = ParseArgs(details)
  ShowVersion(includeArgs=1)

  models = []
  if details.note and details.dbName:
    tblName = extras[0]
    message('-> Retrieving models from database')
    conn = DbConnect(details.dbName, tblName)
    blobs = conn.GetData(fields='model', where="where note='%s'" % (details.note))
    for blob in blobs:
      blob = blob[0]
      try:
        models.append(pickle.loads(str(blob)))
      except Exception:
        import traceback
        traceback.print_exc()
        message('Model load failed')

  else:
    message('-> Loading model')
    modelFile = open(extras[0], 'rb')
    models.append(pickle.load(modelFile))
  if not len(models):
    error('No composite models found')
    sys.exit(-1)
  else:
    message('-> Working with %d models.' % len(models))

  extras = extras[1:]

  for fName in extras:
    if details.dbName != '':
      details.tableName = fName
      data = details.GetDataSet(pickleCol=details.pickleCol,
                                pickleClass=DataStructs.ExplicitBitVect)
    else:
      data = DataUtils.BuildDataSet(fName)
    descNames = data.GetVarNames()
    nModels = len(models)
    screenResults = [None] * nModels
    dataSets = [None] * nModels
    message('-> Constructing and screening data sets')
    testIdx = list(range(data.GetNPts()))
    trainIdx = testIdx

    for modelIdx in range(nModels):
      # tmpD = copy.deepcopy(data)
      tmpD = data
      model = models[modelIdx]
      message('.', noRet=1)

      try:
        seed = model._randomSeed
      except AttributeError:
        pass
      else:
        DataUtils.InitRandomNumbers(seed)

      if details.shuffleActivities or details.randomActivities:
        shuffle = details.shuffleActivities
        randomize = 1
        DataUtils.RandomizeActivities(tmpD, shuffle=details.shuffleActivities, runDetails=details)
      else:
        randomize = False
        shuffle = False

      if hasattr(model, '_shuffleActivities') and \
         model._shuffleActivities and \
         not shuffle:
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
        message('******  WARNING: Shuffled model being screened with unshuffled data.')
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
      if hasattr(model, '_randomizeActivities') and \
         model._randomizeActivities and \
         not randomize:
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')
        message('******  WARNING: Random model being screened with non-random data.')
        message('*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*')

      trainIdx, testIdx = PrepareDataFromDetails(model, details, tmpD, verbose=1)
      screenResults[modelIdx] = CollectResults(testIdx, tmpD, model,
                                               errorEstimate=details.errorEstimate)
      dataSets[modelIdx] = testIdx
    for tol in details.screenVoteTol:
      if len(details.screenVoteTol) > 1:
        message('\n-----*****-----*****-----*****-----*****-----*****-----*****-----\n')
        message('Tolerance: %f' % tol)
      nGood = numpy.zeros(nModels, numpy.float)
      nBad = numpy.zeros(nModels, numpy.float)
      nSkip = numpy.zeros(nModels, numpy.float)
      confGood = numpy.zeros(nModels, numpy.float)
      confBad = numpy.zeros(nModels, numpy.float)
      confSkip = numpy.zeros(nModels, numpy.float)
      if details.enrichTgt >= 0:
        enrichments = numpy.zeros(nModels, numpy.float)
      goodVoteDict = {}
      badVoteDict = {}
      noVoteDict = {}
      voteTab = None
      for modelIdx in range(nModels):
        model = models[modelIdx]
        model.SetInputOrder(descNames)
        testIdx = dataSets[modelIdx]
        screenRes = screenResults[modelIdx]
        if not details.detailedScreen:
          g, b, s, aG, aB, aS, vT = ScreenIt(model, testIdx, tmpD, details.partialVote, tol,
                                             verbose=details.verbose, screenResults=screenRes)
        else:
          if model.GetActivityQuantBounds():
            nRes = len(model.GetActivityQuantBounds()) + 1
          else:
            nRes = model.GetQuantBounds()[1][-1]
          badVotes = []
          noVotes = []
          if (hasattr(details, 'showAll') and details.showAll) or \
             (hasattr(details, 'predPlot') and details.predPlot):
            goodVotes = []
          else:
            goodVotes = None
          g, b, s, aG, aB, aS, vT = ShowVoteResults(
            testIdx, tmpD, model, nRes, tol, verbose=details.verbose, screenResults=screenRes,
            badVotes=badVotes, noVotes=noVotes, goodVotes=goodVotes,
            errorEstimate=details.errorEstimate)
          if voteTab is None:
            voteTab = numpy.zeros(vT.shape, numpy.float)
          if details.errorAnalysis:
            for a, p, c, idx in badVotes:
              label = testIdx[idx]
              if hasattr(details, 'enrichTgt') and details.enrichTgt >= 0:
                if a == details.enrichTgt:
                  badVoteDict[label] = badVoteDict.get(label, 0) + 1
              else:
                badVoteDict[label] = badVoteDict.get(label, 0) + 1
            for a, p, c, idx in noVotes:
              label = testIdx[idx]
              if hasattr(details, 'enrichTgt') and details.enrichTgt >= 0:
                if a == details.enrichTgt:
                  noVoteDict[label] = noVoteDict.get(label, 0) + 1
              else:
                noVoteDict[label] = noVoteDict.get(label, 0) + 1

          if hasattr(details, 'showAll') and details.showAll:
            for a, p, c, idx in goodVotes:
              label = testIdx[idx]
              if details.enrichTgt >= 0:
                if a == details.enrichTgt:
                  goodVoteDict[label] = goodVoteDict.get(label, 0) + 1
              else:
                goodVoteDict[label] = goodVoteDict.get(label, 0) + 1

          if details.enrichTgt > -1:
            enrichments[modelIdx] = CalcEnrichment(vT, tgt=details.enrichTgt)

          voteTab += vT
        if details.detailedScreen and hasattr(details, 'predPlot') and details.predPlot:
          MakePredPlot(details, testIdx, tmpD, goodVotes, badVotes, nRes, verbose=1)

        if hasattr(details, 'showAll') and details.showAll:
          print('-v-v-v-v-v-v-v-    All Votes      -v-v-v-v-v-v-v-')
          print('id, prediction, confidence, flag(-1=skipped,0=wrong,1=correct)')
          for ans, pred, conf, idx in goodVotes:
            pt = tmpD[testIdx[idx]]
            assert model.GetActivityQuantBounds() or pt[-1] == ans, 'bad point?: %s != %s' % (
              str(pt[-1]), str(ans))
            print('%s, %d, %.4f, 1' % (str(pt[0]), pred, conf))
          for ans, pred, conf, idx in badVotes:
            pt = tmpD[testIdx[idx]]
            assert model.GetActivityQuantBounds() or pt[-1] == ans, 'bad point?: %s != %s' % (
              str(pt[-1]), str(ans))
            print('%s, %d, %.4f, 0' % (str(pt[0]), pred, conf))
          for ans, pred, conf, idx in noVotes:
            pt = tmpD[testIdx[idx]]
            assert model.GetActivityQuantBounds() or pt[-1] == ans, 'bad point?: %s != %s' % (
              str(pt[-1]), str(ans))
            print('%s, %d, %.4f, -1' % (str(pt[0]), pred, conf))
          print('-^-^-^-^-^-^-^-  -^-^-^-^-^-^-^-')

        nGood[modelIdx] = g
        nBad[modelIdx] = b
        nSkip[modelIdx] = s
        confGood[modelIdx] = aG
        confBad[modelIdx] = aB
        confSkip[modelIdx] = aS
      print()

      if nModels > 1:
        print('-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')
        print('AVERAGES:')

        avgNBad = sum(nBad) / nModels
        devNBad = numpy.sqrt(sum((nBad - avgNBad)**2) / (nModels - 1))

        bestIdx = numpy.argsort(nBad)[0]

        avgNGood = sum(nGood) / nModels
        devNGood = numpy.sqrt(sum((nGood - avgNGood)**2) / (nModels - 1))

        avgNSkip = sum(nSkip) / nModels
        devNSkip = numpy.sqrt(sum((nSkip - avgNSkip)**2) / (nModels - 1))

        avgConfBad = sum(confBad) / nModels
        devConfBad = numpy.sqrt(sum((confBad - avgConfBad)**2) / (nModels - 1))

        avgConfGood = sum(confGood) / nModels
        devConfGood = numpy.sqrt(sum((confGood - avgConfGood)**2) / (nModels - 1))

        avgConfSkip = sum(confSkip) / nModels
        devConfSkip = numpy.sqrt(sum((confSkip - avgConfSkip)**2) / (nModels - 1))

        nClassified = avgNGood + avgNBad
        nExamples = nClassified + avgNSkip
        print('Misclassifications: \t%%%5.2f(%%%5.2f)   %4.1f(%4.1f) / %d' %
              (100 * avgNBad / nExamples, 100 * devNBad / nExamples, avgNBad, devNBad, nExamples))
        if avgNSkip > 0:
          print('\tthreshold: \t%%%5.2f(%%%5.2f)   %4.1f(%4.1f) / %d' %
                (100 * avgNBad / nClassified, 100 * devNBad / nClassified, avgNBad, devNBad,
                 nClassified))
          print()
          print('Number Skipped: %%%4.2f(%%%4.2f)    %4.2f(%4.2f)' %
                (100 * avgNSkip / nExamples, 100 * devNSkip / nExamples, avgNSkip, devNSkip))

        print()
        print('Confidences:')
        print('\tCorrect: \t%4.2f(%4.2f)' % (100 * avgConfGood, 100 * devConfGood))
        print('\tIncorrect: \t%4.2f(%4.2f)' % (100 * avgConfBad, 100 * devConfBad))
        if avgNSkip > 0:
          print('\tSkipped: \t%4.2f(%4.2f)' % (100 * avgConfSkip, 100 * devConfSkip))

        if details.detailedScreen:
          message('Results Table:')
          voteTab = numpy.transpose(voteTab) / nModels
          nResultCodes = len(voteTab)
          colCounts = numpy.sum(voteTab, 0)
          rowCounts = numpy.sum(voteTab, 1)
          print()
          for i in range(nResultCodes):
            if rowCounts[i] == 0:
              rowCounts[i] = 1
            row = voteTab[i]
            message('    ', noRet=1)
            for j in range(nResultCodes):
              entry = row[j]
              message(' % 6.2f' % entry, noRet=1)
            message('     | % 4.2f' % (100. * voteTab[i, i] / rowCounts[i]))
          message('    ', noRet=1)
          for i in range(nResultCodes):
            message('-------', noRet=1)
          message('')
          message('    ', noRet=1)
          for i in range(nResultCodes):
            if colCounts[i] == 0:
              colCounts[i] = 1
            message(' % 6.2f' % (100. * voteTab[i, i] / colCounts[i]), noRet=1)
          message('')
          if details.enrichTgt > -1:
            mean = sum(enrichments) / nModels
            enrichments -= mean
            dev = numpy.sqrt(sum(enrichments * enrichments)) / (nModels - 1)
            message('   Enrichment of value %d: %.4f (%.4f)' % (details.enrichTgt, mean, dev))
      else:
        bestIdx = 0
      print('------------------------------------------------')
      print('Best Model: ', bestIdx + 1)
      bestBad = nBad[bestIdx]
      bestGood = nGood[bestIdx]
      bestSkip = nSkip[bestIdx]
      nClassified = bestGood + bestBad
      nExamples = nClassified + bestSkip
      print('Misclassifications: \t%%%5.2f   %d / %d' % (100 * bestBad / nExamples, bestBad,
                                                         nExamples))
      if bestSkip > 0:
        print('\tthreshold: \t%%%5.2f   %d / %d' % (100 * bestBad / nClassified, bestBad,
                                                    nClassified))
        print()
        print('Number Skipped: %%%4.2f    %d' % (100 * bestSkip / nExamples, bestSkip))

      print()
      print('Confidences:')
      print('\tCorrect: \t%4.2f' % (100 * confGood[bestIdx]))
      print('\tIncorrect: \t%4.2f' % (100 * confBad[bestIdx]))
      if bestSkip > 0:
        print('\tSkipped: \t%4.2f' % (100 * confSkip[bestIdx]))

      if nModels == 1 and details.detailedScreen:
        message('')
        message('Results Table:')
        voteTab = numpy.transpose(vT)
        nResultCodes = len(vT)
        colCounts = numpy.sum(voteTab, 0)
        rowCounts = numpy.sum(voteTab, 1)
        message('')
        for i in range(nResultCodes):
          if rowCounts[i] == 0:
            rowCounts[i] = 1
          row = voteTab[i]
          message('    ', noRet=1)
          for j in range(nResultCodes):
            entry = row[j]
            message(' % 6.2f' % entry, noRet=1)
          message('     | % 4.2f' % (100. * voteTab[i, i] / rowCounts[i]))
        message('    ', noRet=1)
        for i in range(nResultCodes):
          message('-------', noRet=1)
        message('')
        message('    ', noRet=1)
        for i in range(nResultCodes):
          if colCounts[i] == 0:
            colCounts[i] = 1
          message(' % 6.2f' % (100. * voteTab[i, i] / colCounts[i]), noRet=1)
        message('')
      if details.errorAnalysis:
        message('\n*-*-*-*-*-*-*-*- ERROR ANALYSIS -*-*-*-*-*-*-*-*\n')
        ks = badVoteDict.keys()
        if len(ks):
          message(' ---> Bad Vote Counts')
        ks = noVoteDict.keys()
        if len(ks):
          message(' ---> Skipped Compound Counts')
          for k in ks:
            pt = data[k]
            message('%s,%d' % (str(pt[0]), noVoteDict[k]))

      if hasattr(details, 'showAll') and details.showAll:
        ks = goodVoteDict.keys()
        if len(ks):
          message(' ---> Good Vote Counts')
          for k in ks:
            pt = data[k]
            message('%s,%d' % (str(pt[0]), goodVoteDict[k]))
