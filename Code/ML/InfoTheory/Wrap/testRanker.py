import unittest
import numpy
import os

from rdkit.six.moves import cPickle

from rdkit import RDConfig,RDRandom
from rdkit.ML.InfoTheory import rdInfoTheory as rdit
from rdkit import DataStructs

def feq(a,b,tol=1e-4):
    return abs(a-b)<tol

class TestCase(unittest.TestCase):

    def setUp(self) :
        pass

    def test0GainFuns(self):
       arr = numpy.array([9,5])
       self.assertTrue(feq(rdit.InfoEntropy(arr),0.9403))
       arr = numpy.array([9,9])
       self.assertTrue(feq(rdit.InfoEntropy(arr),1.0000))
       arr = numpy.array([5,5])
       self.assertTrue(feq(rdit.InfoEntropy(arr),1.0000))
       arr = numpy.array([5,0])
       self.assertTrue(feq(rdit.InfoEntropy(arr),0.0000))
       arr = numpy.array([5,5,5])
       self.assertTrue(feq(rdit.InfoEntropy(arr),1.5850))
       arr = numpy.array([2,5,5])
       self.assertTrue(feq(rdit.InfoEntropy(arr),1.4834))

       
       mat2 = numpy.array([[6,2],[3,3]])
       self.assertTrue(feq(rdit.InfoGain(mat2),0.0481))
       self.assertTrue(feq(rdit.ChiSquare(mat2),0.9333))
       
       mat3 = numpy.array([[1,1],[2,1]])
       self.assertTrue(feq(rdit.InfoGain(mat3),0.0200))

       
       mat4 = numpy.array([[2,0],[1,2]])
       self.assertTrue(feq(rdit.InfoGain(mat4),0.4200))


       mat5 = numpy.array([[0,0],[0,0]])
       self.assertTrue(feq(rdit.InfoGain(mat5),0.0000))


       mat6 = numpy.array([[1,0],[1,0]])
       self.assertTrue(feq(rdit.InfoGain(mat6),0.0000))


       
       
    def test1ranker(self) :
        nbits = 100
        ninst = 100
        dm = 50
        nact = 10
        nc = 2
        rn = rdit.InfoBitRanker(nbits, nc, rdit.InfoType.ENTROPY)
        fps = []
        na = 0
        ni = 0
        for i in range(ninst) :
            v = DataStructs.SparseBitVect(nbits)
            for j in range(dm):
                v.SetBit(RDRandom.randrange(0,nbits))

            
            if (RDRandom.randrange(0,ninst) < nact) :
                na += 1
                rn.AccumulateVotes(v, 1)
                fps.append((v,1))
            else:
                ni += 1
                rn.AccumulateVotes(v, 0)
                fps.append((v,0))
                
        res =  rn.GetTopN(50)

        rn2 = rdit.InfoBitRanker(nbits, nc)
        for fp in fps:
            rn2.AccumulateVotes(fp[0], fp[1])

        res2 = rn2.GetTopN(50)
        self.assertTrue((res==res2).all())
        
        rn3 = rdit.InfoBitRanker(nbits, nc, rdit.InfoType.BIASENTROPY)
        #rn3.SetBiasList([0])
        for fp in fps:
            rn3.AccumulateVotes(fp[0], fp[1])

        res3 = rn3.GetTopN(50)
        for i in range(50) :
            fan = res3[i,2]/na
            fin = res3[i,3]/ni
            self.assertTrue(fan > fin)
                          
    def test2ranker(self) :
        nbits = 100
        ninst = 100
        dm = 50
        nact = 10
        nc = 2
        RDRandom.seed(23)
        rn = rdit.InfoBitRanker(nbits, nc, rdit.InfoType.ENTROPY)
        rn.SetMaskBits([63,70,15,25,10])
        fps = []
        na = 0
        ni = 0
        for i in range(ninst) :
            v = DataStructs.SparseBitVect(nbits)
            for j in range(dm):
                v.SetBit(RDRandom.randrange(0,nbits))
            if (RDRandom.randrange(0,ninst) < nact) :
                na += 1
                rn.AccumulateVotes(v, 1)
                fps.append((v,1))
            else:
                ni += 1
                rn.AccumulateVotes(v, 0)
                fps.append((v,0))
        res =  rn.GetTopN(5)
        ids = [int(x[0]) for x in res]
        ids.sort()
        self.assertTrue(ids==[10,15,25,63,70])
        try:
            res = rn.GetTopN(10)
        except:
            ok = 1
        else:
            ok = 0
        self.assertTrue(ok)

    def test3Issue140(self) :
        nbits = 2
        examples = [[0,0,0],[1,1,0],[0,0,1],[1,1,1]]
        rn = rdit.InfoBitRanker(2,2,rdit.InfoType.ENTROPY)
        for example in examples:
            act = example.pop(-1)
            bv = DataStructs.ExplicitBitVect(2)
            for i in range(2):
                bv[i] = example[i]
            rn.AccumulateVotes(bv,act)
        try:
            res =  rn.GetTopN(1)
        except:
            res = None
        self.assertTrue(res is not None)    

    def test4Issue237(self) :
        with open(os.path.join(RDConfig.RDBaseDir,'Code','ML','InfoTheory','Wrap','testData','Issue237.pkl'),'rb') as inF:
            examples,avail,bias,nB,nPoss = cPickle.load(inF, encoding='bytes')
        ranker = rdit.InfoBitRanker(nB,nPoss,rdit.InfoType.BIASENTROPY)
        ranker.SetMaskBits(avail)
        for ex in examples:
            ranker.AccumulateVotes(ex[1],ex[-1])
        # this dumps core on linux if the bug isn't fixed:
        v=ranker.GetTopN(1)
        self.assertTrue(int(v[0][0])==12)
                          
if __name__ == '__main__':
    unittest.main()
                       
