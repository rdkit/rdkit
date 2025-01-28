import numpy

from rdkit.ML.Cluster import Murtagh

print('1')
d = numpy.array([[10.0, 5.0], [20.0, 20.0], [30.0, 10.0], [30.0, 15.0], [5.0, 10.0]], float)
print('2')
# clusters = Murtagh.ClusterData(d,len(d),Murtagh.WARDS)
# for i in range(len(clusters)):
#   clusters[i].Print()
# print('3')

dists = []
for i in range(len(d)):
  for j in range(i):
    dist = sum((d[i] - d[j])**2)
    dists.append(dist)
dists = numpy.array(dists)

print('Wards:')
clusters = Murtagh.ClusterData(dists, len(d), Murtagh.WARDS, isDistData=1)
clusters[0].Print()

print('SLINK:')
clusters = Murtagh.ClusterData(dists, len(d), Murtagh.SLINK, isDistData=1)
clusters[0].Print()

print('CLINK:')
clusters = Murtagh.ClusterData(dists, len(d), Murtagh.CLINK, isDistData=1)
clusters[0].Print()

print('UPGMA:')
clusters = Murtagh.ClusterData(dists, len(d), Murtagh.UPGMA, isDistData=1)
clusters[0].Print()
