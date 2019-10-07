#
#  Copyright (C) 2001-2008  greg Landrum
#
""" unit testing code for the Stats module

"""
import unittest
from rdkit.ML.Data import Stats
import numpy

FLOAT_TOL = 1e-4


def feq(x, y, tol=FLOAT_TOL):
  return abs(x - y) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    self.d = numpy.array([[7, 4, 3], [4, 1, 8], [6, 3, 5], [8, 6, 1], [8, 5, 7], [7, 2, 9],
                          [5, 3, 3], [9, 5, 8], [7, 4, 5], [8, 2, 2]], 'd')

  def testCovariance(self):
    # """ test the covariance x calculation
    #  test case from http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
    # """

    d = numpy.array([[4., 2, 0.6], [4.2, 2.1, 0.59], [3.9, 2.0, 0.58], [4.3, 2.1, 0.62],
                     [4.1, 2.2, 0.63]])
    m = Stats.FormCovarianceMatrix(d)
    target = numpy.array([[0.025, 0.0075, 0.00175], [0.0075, 0.007, 0.00135],
                          [0.00175, 0.00135, 0.00043]])
    diff = abs(m - target)
    assert max(diff.ravel()) < FLOAT_TOL, 'covariance matrix incorrect'

  def testCorrelation(self):
    # " test the correlation matrix calculation "
    m = Stats.FormCorrelationMatrix(self.d)
    target = numpy.array([[1., 0.66865732, -0.10131374], [0.66865732, 1., -0.28792771],
                          [-0.10131374, -0.28792771, 1.]])
    diff = abs(m - target)
    assert max(diff.ravel()) < FLOAT_TOL, 'correlation matrix incorrect'

  def testPCA(self):
    # " test the PCA calculation "
    eVals, eVects = Stats.PrincipalComponents(self.d)
    tVals = numpy.array([1.76877414, 0.92707592, 0.30414995])
    tVect = numpy.array([[-0.64200458, -0.66321742, 0.38467229],
                         [0.34166917, 0.20166619, 0.91792861],
                         [-0.68636164, 0.72074503, 0.09713033]])
    assert max((abs(eVals - tVals)).ravel()) < FLOAT_TOL, 'bad variances from PCA'
    for i in range(eVects.shape[0]):
      assert (max((abs(eVects[i] - tVect[i])).ravel()) < FLOAT_TOL or
              max((abs(eVects[i] + tVect[i])).ravel()) < FLOAT_TOL), 'bad vectors from PCA'

  def testTransform(self):
    # " test transformation to PCA frame "
    _, eVects = Stats.PrincipalComponents(self.d)
    pts = Stats.TransformPoints(eVects, self.d)
    refPs = [numpy.array([-1.20362098, -1.79265006, 0.08776266]),
             numpy.array([4.63540648, 1.1669869, 0.47026415]),
             numpy.array([0.8709456, -0.50012821, 0.24763993]),
             numpy.array([-3.94140499, -2.88350573, 0.64863041]),
             numpy.array([-0.97015382, 2.42239972, 0.51066736]),
             numpy.array([2.43084762, 3.3115892, -0.77094542]),
             numpy.array([0.74360559, -2.67765459, 0.73974091]),
             numpy.array([-1.2274861, 3.6819975, -0.07856395]),
             numpy.array([-0.4342764, 0.04320715, 0.28202332]),
             numpy.array([-0.903863, -2.77224188, -2.13721937])]

    p0 = refPs[0]
    p3 = refPs[3]
    p5 = refPs[5]
    assert max(abs(pts[0] - p0)) < FLOAT_TOL, 'p0 comparison failed %s!=%s' % (str(pts[0]), str(p0))
    assert max(abs(pts[3] - p3)) < FLOAT_TOL, 'p3 comparison failed %s!=%s' % (str(pts[3]), str(p3))
    assert max(abs(pts[5] - p5)) < FLOAT_TOL, 'p5 comparison failed %s!=%s' % (str(pts[5]), str(p5))

  def testTransform2(self):
    # """ testing that rotation of points into PCA frame
    # doesn't change dot products """
    self.d = numpy.array([
      [1, 2.068703704, 2.040555556, 2.068703704, 2.141782407, 7.46],
      [2, 1.48537037, -0.186756425, 1.48537037, 1.803819444, 8.16],
      [3, 1.917469136, 0.785465797, 1.917469136, 2.046875, 8.68],
      [4, 2.068703704, 1.125743575, 2.068703704, 2.131944444, 8.89],
      [5, 2.138703704, 1.283243575, 2.138703704, 2.171319444, 9.25],
      [6, 2.152037037, 1.313243575, 2.152037037, 2.178819444, 9.3],
      [7, 1.730740741, 1.457222222, -0.179901738, 1.558449074, 7.52],
      [8, 1.973796296, 1.889320988, 0.792320484, 1.99054784, 8.16],
      [9, 2.058865741, 2.040555556, 1.132598262, 2.141782407, 8.3],
      [10, 2.098240741, 2.110555556, 1.290098262, 2.211782407, 8.4],
      [11, 2.105740741, 2.123888889, 1.320098262, 2.225115741, 8.46],
      [12, 1.390462963, -0.37502803, 0.171950113, 1.652584877, 8.19],
      [13, 1.475532407, -0.223793462, 0.512227891, 1.803819444, 8.57],
      [14, 1.522407407, -0.140460128, 0.699727891, 1.887152778, 8.82],
      [15, 1.822561728, 0.597194192, 0.604048879, 1.895640432, 8.89],
      [16, 1.907631173, 0.74842876, 0.944326657, 2.046875, 8.92],
      [17, 1.954506173, 0.831762094, 1.131826657, 2.130208333, 8.96],
      [18, 1.973796296, 0.93747197, 0.755283447, 1.980709877, 9],
      [19, 2.058865741, 1.088706538, 1.095561224, 2.131944444, 9.35],
      [20, 2.105740741, 1.172039872, 1.283061224, 2.215277778, 9.22],
      [21, 2.189074074, 1.359539872, 1.366394558, 2.262152778, 9.3],
      [22, 2.142199074, 1.276206538, 1.178894558, 2.178819444, 9.52]
    ], 'd')
    self.d = self.d[:, 1:-2]
    _, eVects = Stats.PrincipalComponents(self.d)
    pts = Stats.TransformPoints(eVects, self.d)

    avg = sum(self.d) / len(self.d)
    self.d -= avg

    for i in range(len(pts)):
      for j in range(len(pts)):
        vi = self.d[i]
        vi /= numpy.sqrt(numpy.dot(vi, vi))
        vj = self.d[j]
        vj /= numpy.sqrt(numpy.dot(vj, vj))
        pvi = pts[i]
        pvi /= numpy.sqrt(numpy.dot(pvi, pvi))
        pvj = pts[j]
        pvj /= numpy.sqrt(numpy.dot(pvj, pvj))
        assert feq(numpy.dot(vi, vj), numpy.dot(pvi, pvj)), 'bad dot: %4.4f %4.4f' % (numpy.dot(
          vi, vj), numpy.dot(pvi, pvj))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
