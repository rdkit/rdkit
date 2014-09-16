from __future__ import print_function
from rdkit import RDConfig
RDConfig.usePgSQL=0
import unittest
from rdkit.ML import InfoTheory
from rdkit import DataStructs
from rdkit.Dbase.DbConnection import DbConnect
import os
from rdkit.six.moves import cPickle as pickle

def feq(v1,v2,tol2=1e-4):
    return abs(v1-v2)<=tol2

def getFingerprints(conn) :
    data = conn.GetData(table='signatures', fields='mol_name,fingerprint')
    fpMap = {}
    for dat in data :
        pkl = str(dat[1])
        sbv = pickle.loads(pkl)
        fpMap[dat[0]] = sbv
    return fpMap

def getNameAct(conn):
    data = conn.GetData(table='raw_data', fields='mol_name,activity_class')
    nameAct = {}
    for dat in data :
        nameAct[dat[0]] = dat[1]

    return nameAct

def ReadCombiInfo(fileName) :
    infil = open(fileName, 'r')
    lines = infil.readlines()
    infil.close()
    infos = []
    for lin in lines:
        tlst = lin.strip().split()
        info = float(tlst[1])
        infos.append(info)
    return infos

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test0Ranker(self) :
        nbits = 5000
        dbName = os.path.join('../','test_data', 'FEW_CDK2.GDB')
        conn = DbConnect(dbName)
        fps = getFingerprints(conn)
        nameAct = getNameAct(conn)
        sl = len(fps.values()[0])
        rnkr = InfoTheory.InfoBitRanker(sl, 2, InfoTheory.InfoType.ENTROPY)
        
        print("Collecting Votes ....")
        for key in nameAct.keys() :
            if nameAct[key] == 100 :
                rnkr.AccumulateVotes(fps[key], 0)
            if nameAct[key] == 0 :
                rnkr.AccumulateVotes(fps[key], 1)

        # now do the ranking
        print("ranking bits ....")
        topN = rnkr.GetTopN(nbits)
        
        # get the combichem ranked list from a file
        cfile = os.path.join('test_data', 'combiRank.out')
        combiInfo = ReadCombiInfo(cfile)
        # now check if the infocontents are the same as the combichem stuff
        print("Comparing bit info contents ....")
        for i in range(900) :
            assert feq(topN[i,1], combiInfo[i])

        ofile = os.path.join('test_data', 'rdTopBits.txt')
        rnkr.WriteTopBitsToFile(ofile)
        
    def test1BiasRanker(self) :
        nbits = 5000
        dbName = os.path.join('../','test_data', 'FEW_CDK2.GDB')
        conn = DbConnect(dbName)
        fps = getFingerprints(conn)
        nameAct = getNameAct(conn)
        sl = len(fps.values()[0])
        rnkr = InfoTheory.InfoBitRanker(sl, 2, InfoTheory.InfoType.BIASENTROPY)
        rnkr.SetBiasList([0])
        print("Collecting Votes ....")
        for key in nameAct.keys() :
            if nameAct[key] == 100 :
                rnkr.AccumulateVotes(fps[key], 0)
            if nameAct[key] == 0 :
                rnkr.AccumulateVotes(fps[key], 1)

        # now do the ranking
        print("ranking bits ....")
        topN = rnkr.GetTopN(nbits)

        # get the combichem ranked list from a file
        cfile = os.path.join('test_data', 'combiRank.out')
        combiInfo = ReadCombiInfo(cfile)
        # now check if the infocontents are the same as the combichem stuff
        print("Comparing bit info contents ....")
        for i in range(nbits) :
            assert feq(topN[i,1], combiInfo[i])
            
    def test2ChiSquare(self) :
        nbits = 5000
        dbName = os.path.join('../','test_data', 'FEW_CDK2.GDB')
        conn = DbConnect(dbName)
        fps = getFingerprints(conn)
        nameAct = getNameAct(conn)
        sl = len(fps.values()[0])
        rnkr = InfoTheory.InfoBitRanker(sl, 2, InfoTheory.InfoType.BIASCHISQUARE)
        rnkr.SetBiasList([0])
        print("Collecting Votes ....")
        for key in nameAct.keys() :
            if nameAct[key] == 100 :
                rnkr.AccumulateVotes(fps[key], 0)
            if nameAct[key] == 0 :
                rnkr.AccumulateVotes(fps[key], 1)

        # now do the ranking
        print("ranking bits ....")
        topN = rnkr.GetTopN(nbits)

        # get the combichem ranked list from a file
        cfile = os.path.join('test_data', 'combiRankChi.out')
        combiInfo = ReadCombiInfo(cfile)
        # now check if the infocontents are the same as the combichem stuff
        print("Comparing bit info contents ....")
        for i in range(nbits) :
            assert feq(topN[i,1], combiInfo[i])
        #rnkr.WriteTopBitsToFile("chiBitsBias.txt")
                
if __name__ == '__main__':
    unittest.main()
    
        
