#
#  Copyright (C) 2000-2008  Greg Landrum and Rational Discovery LLC
#

from rdkit.SimDivFilters import rdSimDivPickers as rdsimdiv
if rdsimdiv is None:
   raise ImportError('rdSimDivPickers not built')
from rdkit import DataStructs
import numpy

class BitClusterer(object):
    """ Class to cluster a set of bits based on their correllation

    The correlation matrix is first built using by reading the fingerprints
    from a database or a list of fingerprints
    """

    def __init__(self, idList, nCluster, type=rdsimdiv.ClusterMethod.WARD):
        self._clusters = []
        self._bidList = idList
        #self._matGen = BitCorrelationMatGenerator(idList)
        self._nClusters = nCluster
        self._type = type

    def ClusterBits(self, corrMat) :
        # clutering code actually needs distances so, take 1/val for each element in corMat
        distMat = 1/corrMat

        pkr = rdsimdiv.HierarchicalClusterPicker(self._type)
        
        cls = pkr.Cluster(distMat, len(self._bidList), self._nClusters)
        # map the clusters to the actual bit ids
        self._clusters = []
        for cl in cls :
            bcls = []
            for i in cl :
                bid = self._bidList[i]
                bcls.append(bid)
            self._clusters.append(bcls)

    def SetClusters(self, clusters):
        assert len(clusters) == self._nClusters
        self._clusters = clusters
        
    def GetClusters(self) :
        return self._clusters

    def MapToClusterScores(self, fp) :
        """ Map the fingerprint to a real valued vector of score based on the bit clusters

        The dimension of the vector is same as the number of clusters. Each value in the 
        vector corresponds to the number of bits in the corresponding cluster
        that are turned on in the fingerprint

        ARGUMENTS:
         - fp : the fingerprint 
        """
    
        scores = [0]*self._nClusters

        i = 0
        for cls in self._clusters:
            for bid in cls :
                if fp[bid] :
                    scores[i] += 1

            i += 1

        return scores

    def MapToClusterFP(self, fp) :
        """ Map the fingerprint to a smaller sized (= number of clusters) fingerprint

        Each cluster get a bit in the new fingerprint and is turned on if any of the bits in
        the cluster are turned on in the original fingerprint"""

        ebv = DataStructs.ExplicitBitVect(self._nClusters)
        i = 0

        for cls in self._clusters:
            for bid in cls :
                if fp[bid] :
                    ebv.SetBit(i)
                    break
            i += 1
        
        return ebv
        
