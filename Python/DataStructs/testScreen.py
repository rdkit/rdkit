## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py

import cPickle,time
from Numeric import *
import numpy.oldnumeric.random_array as RandomArray
from DataStructs import BitRank
nBits = 100
weights = RandomArray.random(nBits)
wSum =  sum(weights)
weights /= wSum
#print 'weights:',weights
nData = 100
thresh = .50
#thresh = .2

dataFile = open('screen.pkl','wb+')

print 'generating data'
d = [None]*nData
acts = zeros(nData,Int)

for i in xrange(nData):
  tmp = RandomArray.randint(0,2,nBits)
  d[i] = tmp
  bitSum = sum(weights*tmp)
  if bitSum > thresh:
    acts[i] = 1

#d = array(d)
#print 'd:'
#print d
#print 'acts:'
#print acts
print 'numActives: %d/%d'%(sum(acts),len(acts))
print 'dumping'
cPickle.dump(d,dataFile)
cPickle.dump(acts,dataFile)

print 'ranking'
t1 = time.time()
bitOrder,infoGains = BitRank.RankBits(d,acts)
t2 = time.time()
print '\t that took %f seconds'%(t2-t1)
#print 'order:',bitOrder
#print 'gains:',infoGains

