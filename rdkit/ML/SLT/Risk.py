#
#  Copyright (C) 2000-2008  greg Landrum
#
""" code for calculating empirical risk

"""
import math


def log2(x):
  return math.log(x) / math.log(2.)


def BurgesRiskBound(VCDim, nData, nWrong, conf):
  """ Calculates Burges's formulation of the risk bound

    The formulation is from Eqn. 3 of Burges's review
    article "A Tutorial on Support Vector Machines for Pattern Recognition"
     In _Data Mining and Knowledge Discovery_ Kluwer Academic Publishers
     (1998) Vol. 2

    **Arguments**

      - VCDim: the VC dimension of the system

      - nData: the number of data points used

      - nWrong: the number of data points misclassified

      - conf: the confidence to be used for this risk bound


    **Returns**

      - a float

    **Notes**

     - This has been validated against the Burges paper

     - I believe that this is only technically valid for binary classification

  """
  #  maintain consistency of notation with Burges's paper
  h = VCDim
  l = nData
  eta = conf

  numerator = h * (math.log(2. * l / h) + 1.) - math.log(eta / 4.)
  structRisk = math.sqrt(numerator / l)

  rEmp = float(nWrong) / l

  return rEmp + structRisk


def CristianiRiskBound(VCDim, nData, nWrong, conf):
  """
    the formulation here is from pg 58, Theorem 4.6 of the book
    "An Introduction to Support Vector Machines" by Cristiani and Shawe-Taylor
    Cambridge University Press, 2000


    **Arguments**

      - VCDim: the VC dimension of the system

      - nData: the number of data points used

      - nWrong: the number of data points misclassified

      - conf: the confidence to be used for this risk bound


    **Returns**

      - a float

    **Notes**

      - this generates odd (mismatching) values

  """
  #  maintain consistency of notation with Christiani's book

  d = VCDim
  delta = conf
  l = nData
  k = nWrong

  structRisk = math.sqrt((4. / nData) * (d * log2((2. * math.e * l) / d) + log2(4. / delta)))
  rEmp = 2. * k / l
  return rEmp + structRisk


def CherkasskyRiskBound(VCDim, nData, nWrong, conf, a1=1.0, a2=2.0):
  """

    The formulation here is from Eqns 4.22 and 4.23 on pg 108 of
    Cherkassky and Mulier's book "Learning From Data" Wiley, 1998.

    **Arguments**

      - VCDim: the VC dimension of the system

      - nData: the number of data points used

      - nWrong: the number of data points misclassified

      - conf: the confidence to be used for this risk bound

      - a1, a2: constants in the risk equation. Restrictions on these values:

          - 0 <= a1 <= 4

          - 0 <= a2 <= 2

    **Returns**

      - a float


    **Notes**

     - This appears to behave reasonably

     - the equality a1=1.0 is by analogy to Burges's paper.

  """
  #  maintain consistency of notation with Cherkassky's book
  h = VCDim
  n = nData
  eta = conf
  rEmp = float(nWrong) / nData

  numerator = h * (math.log(float(a2 * n) / h) + 1) - math.log(eta / 4.)
  eps = a1 * numerator / n

  structRisk = eps / 2. * (1. + math.sqrt(1. + (4. * rEmp / eps)))

  return rEmp + structRisk
