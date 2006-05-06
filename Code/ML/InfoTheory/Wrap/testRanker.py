import RDConfig,RDRandom
import unittest
from ML.InfoTheory import rdInfoTheory as rdit
import DataStructs
from Numeric import *
import os,cPickle

class TestCase(unittest.TestCase):

    def setUp(self) :
        pass

    def test0GainFuns(self):
       arr = array([9,5])
       print arr, rdit.InfoEntropy(arr)
       arr = array([9,9])
       print arr, rdit.InfoEntropy(arr)
       arr = array([5,5])
       print arr, rdit.InfoEntropy(arr)
       arr = array([5,0])
       print arr, rdit.InfoEntropy(arr)
       arr = array([5,5,5])
       print arr, rdit.InfoEntropy(arr)
       arr = array([2,5,5])
       print arr, rdit.InfoEntropy(arr)
       
       mat2 = array([[6,2],[3,3]])
       print 'gain: ', rdit.InfoGain(mat2)
       print "Chi Square: ", rdit.ChiSquare(mat2)
       
       mat3 = array([[1,1],[2,1]])
       print 'gain3: ', rdit.InfoGain(mat3)
       
       mat4 = array([[2,0],[1,2]])
       print 'gain4: ', rdit.InfoGain(mat4)

       mat5 = array([[0,0],[0,0]])
       print 'gain5: ', rdit.InfoGain(mat5)

       mat6 = array([[1,0],[1,0]])
       print 'gain6: ', rdit.InfoGain(mat6)

       
       
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
        assert res == res2
        
        rn3 = rdit.InfoBitRanker(nbits, nc, rdit.InfoType.BIASENTROPY)
        #rn3.SetBiasList([0])
        for fp in fps:
            rn3.AccumulateVotes(fp[0], fp[1])

        res3 = rn3.GetTopN(50)
        for i in range(50) :
            fan = res3[i,2]/na
            fin = res3[i,3]/ni
            assert fan > fin
                          
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
        assert ids==[10,15,25,63,70],str(ids)
        try:
            res = rn.GetTopN(10)
        except:
            ok = 1
        else:
            ok = 0
        assert ok,'expected failure did not happen'

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
        self.failUnless(res is not None)    

    def test4Issue237(self) :
        inF = open(os.path.join(RDConfig.RDBaseDir,'Code','ML','InfoTheory','Wrap','testData','Issue237.pkl'),'rb')
        examples,avail,bias,nB,nPoss = cPickle.load(inF)
        ranker = rdit.InfoBitRanker(nB,nPoss,rdit.InfoType.BIASENTROPY)
        ranker.SetMaskBits(avail)
        for ex in examples:
            ranker.AccumulateVotes(ex[1],ex[-1])
        # this dumps core on linux if the bug isn't fixed:
        v=ranker.GetTopN(1)
        self.failUnless(int(v[0][0])==12)
                          
if __name__ == '__main__':
    unittest.main()
                       
