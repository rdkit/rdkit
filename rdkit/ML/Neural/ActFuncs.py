#
#  Copyright (C) 2000-2008  greg Landrum
#
""" Activation functions for neural network nodes

Activation functions should implement the following API:

 - _Eval(x)_: returns the value of the function at a given point

 - _Deriv(x)_: returns the derivative of the function at a given point

The current Backprop implementation also requires:

 - _DerivFromVal(val)_: returns the derivative of the function when its
                        value is val

In all cases _x_ is a float as is the value returned.

"""
import math


class ActFunc(object):
  """ "virtual base class" for activation functions

  """

  def __call__(self, x):
    return self.Eval(x)


class Sigmoid(ActFunc):
  """ the standard sigmoidal function """

  def Eval(self, x):
    return 1. / (1. + math.exp(-self.beta * x))

  def Deriv(self, x):
    val = self.Eval(x)
    return self.beta * val * (1. - val)

  def DerivFromVal(self, val):
    return self.beta * val * (1. - val)

  def __init__(self, beta=1.):
    self.beta = beta


class TanH(ActFunc):
  """ the standard hyperbolic tangent function """

  def Eval(self, x):
    v1 = math.exp(self.beta * x)
    v2 = math.exp(-self.beta * x)
    return (v1 - v2) / (v1 + v2)

  def Deriv(self, x):
    val = self.Eval(x)
    return self.beta * (1 - val * val)

  def DerivFromVal(self, val):
    return self.beta * (1 - val * val)

  def __init__(self, beta=1.):
    self.beta = beta
