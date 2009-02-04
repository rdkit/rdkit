#
#  Copyright (C) 2000-2008  greg Landrum
#
""" Activation functions for neural network nodes

Activation functions should implement the following API:

 - _Eval(input)_: returns the value of the function at a given point

 - _Deriv(input)_: returns the derivative of the function at a given point

The current Backprop implementation also requires:

 - _DerivFromVal(val)_: returns the derivative of the function when its
                        value is val

In all cases _input_ is a float as is the value returned.

"""
import math
  
class ActFunc(object):
  """ "virtual base class" for activation functions

  """
  def __call__(self,input):
    return self.Eval(input)


class Sigmoid(ActFunc):
  """ the standard sigmoidal function """
  def Eval(self,input):
    return 1./(1.+math.exp(-self.beta*input))

  def Deriv(self,input):
    val = self.Eval(input)
    return self.beta * val * (1. - val)

  def DerivFromVal(self,val):
    return self.beta * val * (1. - val)

  def __init__(self,beta=1.):
    self.beta=beta
    
class TanH(ActFunc):
  """ the standard hyperbolic tangent function """
  def Eval(self,input):
    v1 = math.exp(self.beta*input)
    v2 = math.exp(-self.beta*input)
    return (v1 - v2)/(v1 + v2)
  
  def Deriv(self,input):
    val = self.Eval(input)
    return self.beta * (1 - val*val)

  def DerivFromVal(self,val):
    return self.beta * (1 - val*val)

  def __init__(self,beta=1.):
    self.beta = beta
