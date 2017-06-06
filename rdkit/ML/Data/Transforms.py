#
#  Copyright (C) 2002  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#


def _CenterTForm(dataSet):
  """ INTERNAL USE ONLY

  """
  pass


def _NormalizeTForm(dataSet):
  """ INTERNAL USE ONLY

  """
  pass


def _StandardTForm(dataSet):
  """ INTERNAL USE ONLY

  """
  pass


_availTransforms = [
  ('Center', _CenterTForm, 'translates so that mean(x)=0'),
  ('Normalize', _NormalizeTForm, 'scales so that dot(x,x)=1'),
  ('Standardize', _StandardTForm, 'scales so that dev(x)=0'),
]


def GetAvailTransforms():
  """ returns the list of available data transformations

   **Returns**

     a list of 3-tuples

       1) name of the transform (text)

       2) function describing the transform (should take an
          _MLDataSet_ as an argument)

       3) description of the transform (text)

  """
  global _availTransforms
  return _availTransforms
