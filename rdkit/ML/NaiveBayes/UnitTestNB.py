# $Id$
#
#  Copyright (C) 2004-2005  Rational Discovery LLC
#   All Rights Reserved
#
import sys
import unittest
from rdkit import RDConfig
import os
from rdkit.ML.Data import DataUtils, MLData
from rdkit.ML.NaiveBayes import CrossValidate
from rdkit.DataStructs import ExplicitBitVect


def feq(a,b,tol=1e-4):
    return abs(a-b)<tol

class TestCase(unittest.TestCase):
    def setUp(self) :
        DataUtils.InitRandomNumbers((25,25))

    def test1NaiveBayes(self) :
        fName = os.path.join(RDConfig.RDCodeDir,'ML','NaiveBayes','test_data','stddata.csv')
        data = DataUtils.TextFileToData(fName)
        
        examples = data.GetNamedData()

        nvars = data.GetNVars()
        attrs = range(1,nvars+1)
        npvals= [0] + [3]*nvars + [2]
        qBounds = [0] + [2]*nvars + [0]
        mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds)
        
        assert feq(mod._classProbs[0], 0.54167)
        assert feq(mod._classProbs[1], 0.45833)
        assert feq(mod._QBoundVals[1][0], -0.56995)
        assert feq(mod._QBoundVals[1][1], 0.114)
        assert feq(mod._QBoundVals[2][0], -0.7022)
        assert feq(mod._QBoundVals[2][1], -0.2347)
        assert feq(mod._QBoundVals[3][0], -0.3659)
        assert feq(mod._QBoundVals[3][1], 1.17275)
        
        assert feq(err, 0.16129)

    def test2NaiveBayes(self) :
        fName = os.path.join(RDConfig.RDCodeDir,'ML','NaiveBayes','test_data','stddata.csv')
        data = DataUtils.TextFileToData(fName)
        examples = data.GetNamedData()

        nvars = data.GetNVars()
        attrs = range(1,nvars+1)
        npvals= [0] + [3]*nvars + [2]
        qBounds = [0] + [2]*nvars + [0]
        mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds,
                                                       mEstimateVal=20.0)
        
        assert feq(err, 0.19354)
        
    def test3(self) :
        examples = [
            ['a',  1,0,1,0,  1],
            ['b',  1,0,0,0,  1],
            ['c',  1,0,1,0,  0],
            ['d',  0,1,1,0,  0],
            ['e',  0,1,1,1,  0],
            ]

        nvars = len(examples[0])-2
        attrs = range(1,nvars+1)
        npvals = [0] + [2]*nvars + [2]
        qBounds = [0] + [0]*nvars + [0]
        mdl = CrossValidate.makeNBClassificationModel(examples,attrs,npvals,qBounds)
        nWrong = 0
        for eg in examples:
            p = mdl.ClassifyExample(eg)
            if p != eg[-1]: nWrong +=1
        self.failUnless(nWrong==1)

        bitEx = []
        for eg in examples:
            newEg = [eg[0],None,eg[-1]]
            bv = ExplicitBitVect(nvars)
            for i in range(nvars):
                if eg[i+1]: bv.SetBit(i)
            newEg[1] = bv
            bitEx.append(newEg)
            
        attrs = range(nvars)
        mdl2 = CrossValidate.makeNBClassificationModel(bitEx,attrs,npvals,qBounds,
                                                       useSigs=True)
        nWrong = 0
        for eg in bitEx:
            p = mdl2.ClassifyExample(eg)
            if p != eg[-1]: nWrong +=1
        self.failUnless(nWrong==1)
        
        # now compare:
        for i in range(len(bitEx)):
            eg = examples[i]
            p1 = mdl.ClassifyExample(eg)
            bitEg = bitEx[i]
            p2 = mdl2.ClassifyExample(bitEg)
            self.failUnless(p1==p2)
            v1 = mdl.GetClassificationDetails()
            v2 = mdl.GetClassificationDetails()
            self.failUnless(feq(p1,p2))

    def test4(self) :
        examples = [
            ['a',  1,0,1,0,  1],
            ['b',  1,0,0,0,  1],
            ['c',  1,0,1,0,  0],
            ['d',  0,1,1,0,  0],
            ['e',  0,1,1,1,  0],
            ]

        nvars = len(examples[0])-2
        attrs = range(1,nvars+1)
        origNVars=nvars
        nvars = 10
        npvals = [0] + [2]*nvars + [2]
        qBounds = [0] + [0]*nvars + [0]

        bitEx = []
        for eg in examples:
            newEg = [eg[0],None,eg[-1]]
            bv = ExplicitBitVect(nvars)
            for i in range(origNVars):
                if eg[i+1]: bv.SetBit(i)

            # this bit will yield perfect accuracy if
            #  the attrs argument isn't being used properly:
            if eg[-1]: bv.SetBit(origNVars)
            newEg[1] = bv
            bitEx.append(newEg)
            
        attrs = range(origNVars)
        mdl2 = CrossValidate.makeNBClassificationModel(bitEx,attrs,npvals,qBounds,
                                                       useSigs=True)
        nWrong = 0
        for eg in bitEx:
            p = mdl2.ClassifyExample(eg)
            if p != eg[-1]: nWrong +=1
        self.failUnless(nWrong==1)

    def test5(self) :
        examples = [
            ['a',  1,0,1,0,1,1,0,  1],
            ['b',  1,0,0,0,1,0,0,  1],
            ['c',  1,0,1,0,1,1,0,  0],
            ['d',  0,1,1,0,1,0,0,  0],
            ['e',  0,1,1,1,0,1,0,  0],
            ]

        nvars = len(examples[0])-2
        attrs = range(1,nvars+1)
        npvals = [0] + [2]*nvars + [2]
        qBounds = [0] + [0]*nvars + [0]

        bitEx = []
        for eg in examples:
            newEg = [eg[0],None,eg[-1]]
            bv = ExplicitBitVect(nvars)
            for i in range(nvars):
                if eg[i+1]: bv.SetBit(i)

            # this bit will yield perfect accuracy if
            #  the attrs argument isn't being used properly:
            newEg[1] = bv
            bitEx.append(newEg)
            
        attrs = range(nvars)
        mdl2 = CrossValidate.makeNBClassificationModel(bitEx,attrs,npvals,qBounds,
                                                       useSigs=True,useCMIM=2)
        nWrong = 0
        for eg in bitEx:
            p = mdl2.ClassifyExample(eg)
            if p != eg[-1]: nWrong +=1
        self.failUnless(nWrong==1)
            
if __name__ == '__main__':
    unittest.main()                                                  
  
