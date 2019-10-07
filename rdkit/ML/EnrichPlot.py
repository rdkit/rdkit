# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Command line tool to construct an enrichment plot from saved composite models

Usage:  EnrichPlot [optional args] -d dbname -t tablename <models>

Required Arguments:
  -d "dbName": the name of the database for screening

  -t "tablename": provide the name of the table with the data to be screened

  <models>: file name(s) of pickled composite model(s).
     If the -p argument is also provided (see below), this argument is ignored.

Optional Arguments:
  - -a "list": the list of result codes to be considered active.  This will be
        eval'ed, so be sure that it evaluates as a list or sequence of
        integers. For example, -a "[1,2]" will consider activity values 1 and 2
        to be active

  - --enrich "list": identical to the -a argument above.

  - --thresh: sets a threshold for the plot.  If the confidence falls below
          this value, picking will be terminated

  - -H: screen only the hold out set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -T: screen only the training set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

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

  - -p "tableName": provides the name of a db table containing the
      models to be screened.  If you use this argument, you should also
      use the -N argument (below) to specify a note value.

  - -N "note": provides a note to be used to pull models from a db table.

  - --plotFile "filename": writes the data to an output text file (filename.dat)
    and creates a gnuplot input file (filename.gnu) to plot it

  - --showPlot: causes the gnuplot plot constructed using --plotFile to be
    displayed in gnuplot.

"""
# from rdkit.Dbase.DbConnection import DbConnect


import sys

import numpy

from rdkit import DataStructs
from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import CompositeRun
from rdkit.ML.Data import DataUtils, SplitData, Stats
import pickle


__VERSION_STRING = "2.4.0"


def cmp(t1, t2):
    return (t1 < t2) * -1 or (t1 > t2) * 1


def message(msg, noRet=0, dest=sys.stderr):
    """ emits messages to _sys.stderr_
      override this in modules which import this one to redirect output

      **Arguments**

        - msg: the string to be displayed

    """
    if noRet:
        dest.write('%s ' % (msg))
    else:
        dest.write('%s\n' % (msg))


def error(msg, dest=sys.stderr):
    """ emits messages to _sys.stderr_
      override this in modules which import this one to redirect output

      **Arguments**

        - msg: the string to be displayed

    """
    sys.stderr.write('ERROR: %s\n' % (msg))


def ScreenModel(mdl, descs, data, picking=[1], indices=[], errorEstimate=0):
    """ collects the results of screening an individual composite model that match
      a particular value

     **Arguments**

       - mdl: the composite model

       - descs: a list of descriptor names corresponding to the data set

       - data: the data set, a list of points to be screened.

       - picking: (Optional) a list of values that are to be collected.
         For examples, if you want an enrichment plot for picking the values
         1 and 2, you'd having picking=[1,2].

      **Returns**

        a list of 4-tuples containing:

           - the id of the point

           - the true result (from the data set)

           - the predicted result

           - the confidence value for the prediction

    """
    mdl.SetInputOrder(descs)

    for j in range(len(mdl)):
        tmp = mdl.GetModel(j)
        if hasattr(tmp, '_trainIndices') and not isinstance(tmp._trainIndices, dict):
            tis = {}
            if hasattr(tmp, '_trainIndices'):
                for v in tmp._trainIndices:
                    tis[v] = 1
            tmp._trainIndices = tis

    res = []
    if mdl.GetQuantBounds():
        needsQuant = 1
    else:
        needsQuant = 0

    if not indices:
        indices = list(range(len(data)))
    nTrueActives = 0
    for i in indices:
        if errorEstimate:
            use = []
            for j in range(len(mdl)):
                tmp = mdl.GetModel(j)
                if not tmp._trainIndices.get(i, 0):
                    use.append(j)
        else:
            use = None
        pt = data[i]
        pred, conf = mdl.ClassifyExample(pt, onlyModels=use)
        if needsQuant:
            pt = mdl.QuantizeActivity(pt[:])
        trueRes = pt[-1]
        if trueRes in picking:
            nTrueActives += 1
        if pred in picking:
            res.append((pt[0], trueRes, pred, conf))
    return nTrueActives, res


def AccumulateCounts(predictions, thresh=0, sortIt=1):
    """  Accumulates the data for the enrichment plot for a single model

      **Arguments**

        - predictions: a list of 3-tuples (as returned by _ScreenModels_)

        - thresh: a threshold for the confidence level.  Anything below
          this threshold will not be considered

        - sortIt: toggles sorting on confidence levels


      **Returns**

        - a list of 3-tuples:

          - the id of the active picked here

          - num actives found so far

          - number of picks made so far

    """
    if sortIt:
        predictions.sort(lambda x, y: cmp(y[3], x[3]))
    res = []
    nCorrect = 0
    nPts = 0
    for i in range(len(predictions)):
        ID, real, pred, conf = predictions[i]
        if conf > thresh:
            if pred == real:
                nCorrect += 1
            nPts += 1
            res.append((ID, nCorrect, nPts))

    return res


def MakePlot(details, final, counts, pickVects, nModels, nTrueActs=-1):
    if not hasattr(details, 'plotFile') or not details.plotFile:
        return

    dataFileName = '%s.dat' % (details.plotFile)
    outF = open(dataFileName, 'w+')
    i = 0
    while i < len(final) and counts[i] != 0:
        if nModels > 1:
            _, sd = Stats.MeanAndDev(pickVects[i])
            confInterval = Stats.GetConfidenceInterval(sd, len(pickVects[i]), level=90)
            outF.write('%d %f %f %d %f\n' % (i + 1, final[i][0] / counts[i], final[i][1] / counts[i],
                                             counts[i], confInterval))
        else:
            outF.write('%d %f %f %d\n' % (i + 1, final[i][0] / counts[i], final[i][1] / counts[i],
                                          counts[i]))
        i += 1
    outF.close()
    plotFileName = '%s.gnu' % (details.plotFile)
    gnuF = open(plotFileName, 'w+')
    gnuHdr = """# Generated by EnrichPlot.py version: %s
  set size square 0.7
  set xr [0:]
  set data styl points
  set ylab 'Num Correct Picks'
  set xlab 'Num Picks'
  set grid
  set nokey
  set term postscript enh color solid "Helvetica" 16
  set term X
  """ % (__VERSION_STRING)
    print(gnuHdr, file=gnuF)
    if nTrueActs > 0:
        print('set yr [0:%d]' % nTrueActs, file=gnuF)
    print('plot x with lines', file=gnuF)
    if nModels > 1:
        everyGap = i / 20
        print('replot "%s" using 1:2 with lines,' % (dataFileName), end='', file=gnuF)
        print('"%s" every %d using 1:2:5 with yerrorbars' % (dataFileName, everyGap), file=gnuF)
    else:
        print('replot "%s" with points' % (dataFileName), file=gnuF)
    gnuF.close()

    if hasattr(details, 'showPlot') and details.showPlot:
        try:
            from Gnuplot import Gnuplot
            p = Gnuplot()
            p('load "%s"' % (plotFileName))
            input('press return to continue...\n')
        except Exception:
            import traceback
            traceback.print_exc()


def Usage():
    """ displays a usage message and exits """
    sys.stderr.write(__doc__)
    sys.exit(-1)


if __name__ == '__main__':
    import getopt
    try:
        args, extras = getopt.getopt(sys.argv[1:], 'd:t:a:N:p:cSTHF:v:',
                                     ('thresh=', 'plotFile=', 'showPlot', 'pickleCol=', 'OOB', 'noSort',
                                      'pickBase=', 'doROC', 'rocThresh=', 'enrich='))
    except Exception:
        import traceback
        traceback.print_exc()
        Usage()

    details = CompositeRun.CompositeRun()
    CompositeRun.SetDefaults(details)

    details.activeTgt = [1]
    details.doTraining = 0
    details.doHoldout = 0
    details.dbTableName = ''
    details.plotFile = ''
    details.showPlot = 0
    details.pickleCol = -1
    details.errorEstimate = 0
    details.sortIt = 1
    details.pickBase = ''
    details.doROC = 0
    details.rocThresh = -1
    for arg, val in args:
        if arg == '-d':
            details.dbName = val
        if arg == '-t':
            details.dbTableName = val
        elif arg == '-a' or arg == '--enrich':
            details.activeTgt = eval(val)
            if not isinstance(details.activeTgt, (tuple, list)):
                # if (type(details.activeTgt) not in (types.TupleType, types.ListType)):
                details.activeTgt = (details.activeTgt, )

        elif arg == '--thresh':
            details.threshold = float(val)
        elif arg == '-N':
            details.note = val
        elif arg == '-p':
            details.persistTblName = val
        elif arg == '-S':
            details.shuffleActivities = 1
        elif arg == '-H':
            details.doTraining = 0
            details.doHoldout = 1
        elif arg == '-T':
            details.doTraining = 1
            details.doHoldout = 0
        elif arg == '-F':
            details.filterFrac = float(val)
        elif arg == '-v':
            details.filterVal = float(val)
        elif arg == '--plotFile':
            details.plotFile = val
        elif arg == '--showPlot':
            details.showPlot = 1
        elif arg == '--pickleCol':
            details.pickleCol = int(val) - 1
        elif arg == '--OOB':
            details.errorEstimate = 1
        elif arg == '--noSort':
            details.sortIt = 0
        elif arg == '--doROC':
            details.doROC = 1
        elif arg == '--rocThresh':
            details.rocThresh = int(val)
        elif arg == '--pickBase':
            details.pickBase = val

    if not details.dbName or not details.dbTableName:
        Usage()
        print('*******Please provide both the -d and -t arguments')

    message('Building Data set\n')
    dataSet = DataUtils.DBToData(details.dbName, details.dbTableName, user=RDConfig.defaultDBUser,
                                 password=RDConfig.defaultDBPassword, pickleCol=details.pickleCol,
                                 pickleClass=DataStructs.ExplicitBitVect)

    descs = dataSet.GetVarNames()
    nPts = dataSet.GetNPts()
    message('npts: %d\n' % (nPts))
    final = numpy.zeros((nPts, 2), numpy.float)
    counts = numpy.zeros(nPts, numpy.integer)
    selPts = [None] * nPts

    models = []
    if details.persistTblName:
        conn = DbConnect(details.dbName, details.persistTblName)
        message('-> Retrieving models from database')
        curs = conn.GetCursor()
        curs.execute("select model from %s where note='%s'" %
                     (details.persistTblName, details.note))
        message('-> Reconstructing models')
        try:
            blob = curs.fetchone()
        except Exception:
            blob = None
        while blob:
            message(' Building model %d' % len(models))
            blob = blob[0]
            try:
                models.append(pickle.loads(str(blob)))
            except Exception:
                import traceback
                traceback.print_exc()
                print('Model failed')
            else:
                message('  <-Done')
            try:
                blob = curs.fetchone()
            except Exception:
                blob = None
        curs = None
    else:
        for modelName in extras:
            try:
                model = pickle.load(open(modelName, 'rb'))
            except Exception:
                import traceback
                print('problems with model %s:' % modelName)
                traceback.print_exc()
            else:
                models.append(model)
    nModels = len(models)
    pickVects = {}
    halfwayPts = [1e8] * len(models)
    for whichModel, model in enumerate(models):
        tmpD = dataSet
        try:
            seed = model._randomSeed
        except AttributeError:
            pass
        else:
            DataUtils.InitRandomNumbers(seed)
        if details.shuffleActivities:
            DataUtils.RandomizeActivities(tmpD, shuffle=1)
        if hasattr(model, '_splitFrac') and (details.doHoldout or details.doTraining):
            trainIdx, testIdx = SplitData.SplitIndices(tmpD.GetNPts(), model._splitFrac, silent=1)
            if details.filterFrac != 0.0:
                trainFilt, temp = DataUtils.FilterData(tmpD, details.filterVal, details.filterFrac, -1,
                                                       indicesToUse=trainIdx, indicesOnly=1)
                testIdx += temp
                trainIdx = trainFilt
            if details.doTraining:
                testIdx, trainIdx = trainIdx, testIdx
        else:
            testIdx = list(range(tmpD.GetNPts()))

        message('screening %d examples' % (len(testIdx)))
        nTrueActives, screenRes = ScreenModel(model, descs, tmpD, picking=details.activeTgt,
                                              indices=testIdx, errorEstimate=details.errorEstimate)
        message('accumulating')
        runningCounts = AccumulateCounts(screenRes, sortIt=details.sortIt, thresh=details.threshold)
        if details.pickBase:
            pickFile = open('%s.%d.picks' % (details.pickBase, whichModel + 1), 'w+')
        else:
            pickFile = None

        for i, entry in enumerate(runningCounts):
            entry = runningCounts[i]
            selPts[i] = entry[0]
            final[i][0] += entry[1]
            final[i][1] += entry[2]
            v = pickVects.get(i, [])
            v.append(entry[1])
            pickVects[i] = v
            counts[i] += 1
            if pickFile:
                pickFile.write('%s\n' % (entry[0]))
            if entry[1] >= nTrueActives / 2 and entry[2] < halfwayPts[whichModel]:
                halfwayPts[whichModel] = entry[2]
        message('Halfway point: %d\n' % halfwayPts[whichModel])

    if details.plotFile:
        MakePlot(details, final, counts, pickVects, nModels, nTrueActs=nTrueActives)
    else:
        if nModels > 1:
            print('#Index\tAvg_num_correct\tConf90Pct\tAvg_num_picked\tNum_picks\tlast_selection')
        else:
            print('#Index\tAvg_num_correct\tAvg_num_picked\tNum_picks\tlast_selection')

        i = 0
        while i < nPts and counts[i] != 0:
            if nModels > 1:
                mean, sd = Stats.MeanAndDev(pickVects[i])
                confInterval = Stats.GetConfidenceInterval(sd, len(pickVects[i]), level=90)
                print('%d\t%f\t%f\t%f\t%d\t%s' % (i + 1, final[i][0] / counts[i], confInterval,
                                                  final[i][1] / counts[i], counts[i], str(selPts[i])))
            else:
                print('%d\t%f\t%f\t%d\t%s' % (i + 1, final[i][0] / counts[i], final[i][1] / counts[i],
                                              counts[i], str(selPts[i])))
            i += 1

    mean, sd = Stats.MeanAndDev(halfwayPts)
    print('Halfway point: %.2f(%.2f)' % (mean, sd))
