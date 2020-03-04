# $Id$
#
#  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
#      All Rights Reserved
#
""" Defines Naive Baysean classification model
   Based on development in: Chapter 6 of "Machine Learning" by Tom Mitchell

"""
import numpy
from rdkit.ML.Data import Quantize


def _getBinId(val, qBounds):
    bid = 0
    for bnd in qBounds:
        if (val > bnd):
            bid += 1
    return bid


# FIX: this class has not been updated to new-style classes
# (RD Issue380) because that would break all of our legacy pickled
# data. Until a solution is found for this breakage, an update is
# impossible.
class NaiveBayesClassifier:
    """
    _NaiveBayesClassifier_s can save the following pieces of internal state, accessible via
    standard setter/getter functions:

    1) _Examples_: a list of examples which have been predicted

    2) _TrainingExamples_: List of training examples - the descriptor value of these examples
      are quantized based on info gain using ML/Data/Quantize.py if necessary

    3) _TestExamples_: the list of examples used to test the model

    4) _BadExamples_ : list of examples that were incorrectly classified

    4) _QBoundVals_: Quant bound values for each varaible - a list of lists

    5) _QBounds_ : Number of bounds for each variable

    """

    def __init__(self, attrs, nPossibleVals, nQuantBounds, mEstimateVal=-1.0, useSigs=False):
        """ Constructor

        """
        self._attrs = attrs
        self._mEstimateVal = mEstimateVal
        self._useSigs = useSigs

        self._classProbs = {}

        self._examples = []
        self._trainingExamples = []
        self._testExamples = []
        self._badExamples = []
        self._QBoundVals = {}
        self._nClasses = nPossibleVals[-1]
        self._qBounds = nQuantBounds
        self._nPosVals = nPossibleVals
        self._needsQuant = 1

        self._name = ""
        self.mprob = -1.0

        # for the sake a of efficiency lets try to change the conditional probabilities
        # to a numpy array instead of a dictionary. The three dimension array is indexed
        # on the activity class, the descriptor ID and the descriptor binID
        # self._condProbs = {}
        # self._condProbs = numpy.zeros((self._nClasses, max(self._attrs)+1,
        #                                max(self._nPosVals)+1), 'd')
        self._condProbs = [None] * self._nClasses
        for i in range(self._nClasses):
            if not (hasattr(self, '_useSigs') and self._useSigs):
                nA = max(self._attrs) + 1
                self._condProbs[i] = [None] * nA
                for j in range(nA):
                    nV = self._nPosVals[j]
                    if self._qBounds[j]:
                        nV = max(nV, self._qBounds[j] + 1)
                    self._condProbs[i][j] = [0.0] * nV
            else:
                self._condProbs[i] = {}
                for idx in self._attrs:
                    self._condProbs[i][idx] = [0.0] * 2

    def GetName(self):
        return self._name

    def SetName(self, name):
        self._name = name

    def NameModel(self, varNames):
        self.SetName('NaiveBayesClassifier')

    def GetExamples(self):
        return self._examples

    def SetExamples(self, examples):
        self._examples = examples

    def GetTrainingExamples(self):
        return self._trainingExamples

    def SetTrainingExamples(self, examples):
        self._trainingExamples = examples

    def GetTestExamples(self):
        return self._testExamples

    def SetTestExamples(self, examples):
        self._testExamples = examples

    def SetBadExamples(self, examples):
        self._badExamples = examples

    def GetBadExamples(self):
        return self._badExamples

    def _computeQuantBounds(self):
        neg = len(self._trainingExamples)
        natr = len(self._attrs)

        # make a list of results and values
        allVals = numpy.zeros((neg, natr), 'd')
        res = []  # list of y values
        i = 0
        for eg in self._trainingExamples:
            res.append(eg[-1])
            j = 0
            for ai in self._attrs:
                val = eg[ai]
                allVals[i, j] = val
                j += 1
            i += 1

        # now loop over each of the columns and compute the bounds
        # the number of bounds is determined by the maximum info gain
        i = 0
        for ai in self._attrs:
            nbnds = self._qBounds[ai]
            if nbnds > 0:
                mbnds = []
                mgain = -1.0

                for j in range(1, nbnds + 1):
                    bnds, igain = Quantize.FindVarMultQuantBounds(
                        allVals[:, i], j, res, self._nClasses)
                    if (igain > mgain):
                        mbnds = bnds
                        mgain = igain
                self._QBoundVals[ai] = mbnds
            i += 1

    def trainModel(self):
        """ We will assume at this point that the training examples have been set

        We have to estmate the conditional probabilities for each of the (binned) descriptor
        component give a outcome (or class). Also the probabilities for each class is estimated
        """
        # first estimate the class probabilities
        n = len(self._trainingExamples)
        for i in range(self._nClasses):
            self._classProbs[i] = 0.0

        # for i in range(self._nClasses):
        #   self._classProbs[i] = float(self._classProbs[i])/n

        # first find the bounds for each descriptor value if necessary
        if not self._useSigs and max(self._qBounds) > 0:
            self._computeQuantBounds()

        # now compute the probabilities
        ncls = {}

        incr = 1.0 / n
        for eg in self._trainingExamples:
            cls = eg[-1]
            self._classProbs[cls] += incr
            ncls[cls] = ncls.get(cls, 0) + 1
            tmp = self._condProbs[cls]
            if not self._useSigs:
                for ai in self._attrs:
                    bid = eg[ai]
                    if self._qBounds[ai] > 0:
                        bid = _getBinId(bid, self._QBoundVals[ai])
                    tmp[ai][bid] += 1.0
            else:
                for ai in self._attrs:
                    if eg[1].GetBit(ai):
                        tmp[ai][1] += 1.0
                    else:
                        tmp[ai][0] += 1.0

        # for key in self._condProbs:
        for cls in range(self._nClasses):
            if cls not in ncls:
                continue
            # cls = key[0]
            tmp = self._condProbs[cls]
            for ai in self._attrs:
                if not self._useSigs:
                    nbnds = self._nPosVals[ai]
                    if (self._qBounds[ai] > 0):
                        nbnds = self._qBounds[ai]
                else:
                    nbnds = 2
                for bid in range(nbnds):
                    if self._mEstimateVal <= 0.0:
                        # this is simple the fraction of of time this descriptor component assume
                        # this value for the examples that belong a specific class
                        # self._condProbs[key] = (float(self._condProbs[key]))/ncls[cls]
                        tmp[ai][bid] /= ncls[cls]
                    else:
                        # this a bit more complicated form - more appropriate for unbalanced data
                        # see "Machine Learning" by Tom Mitchell section 6.9.1.1

                        # this is the probability that this descriptor component can take this specific value
                        # in the lack of any other information is is simply the inverse of the number of
                        # possible values 'npossible'
                        # If we quantized this component then
                        #   npossible = 1 + len(self._QBoundVals[ai])
                        # else if we did no qunatize (the descriptor came quantized)
                        #   npossible = nPossibleVals[ai]
                        # ai = key[1]
                        pdesc = 0.0
                        if self._qBounds[ai] > 0:
                            pdesc = 1.0 / (1 + len(self._QBoundVals[ai]))
                        elif (self._nPosVals[ai] > 0):
                            pdesc = 1.0 / (self._nPosVals[ai])
                        else:
                            raise ValueError(
                                'Neither Bounds set nor data pre-quantized for attribute ' + str(ai))
                        tmp[ai][bid] += (self._mEstimateVal) * pdesc
                        tmp[ai][bid] /= (ncls[cls] + self._mEstimateVal)

    def ClassifyExamples(self, examples, appendExamples=0):
        preds = []
        for eg in examples:
            pred = self.ClassifyExample(eg, appendExamples)
            preds.append(int(pred))
        return preds

    def GetClassificationDetails(self):
        """ returns the probability of the last prediction """
        return self.mprob

    def ClassifyExample(self, example, appendExamples=0):
        """ Classify an example by summing over the conditional probabilities
        The most likely class is the one with the largest probability
        """
        if appendExamples:
            self._examples.append(example)
        clsProb = {}
        for key, prob in self._classProbs.items():
            clsProb[key] = prob
            tmp = self._condProbs[key]
            for ai in self._attrs:
                if not (hasattr(self, '_useSigs') and self._useSigs):
                    bid = example[ai]
                    if self._qBounds[ai] > 0:
                        bid = _getBinId(bid, self._QBoundVals[ai])
                else:
                    if example[1].GetBit(ai):
                        bid = 1
                    else:
                        bid = 0
                clsProb[key] *= tmp[ai][bid]

        mkey = -1
        self.mprob = -1.0
        for key, prob in clsProb.items():
            if (prob > self.mprob):
                mkey = key
                self.mprob = prob

        return mkey
