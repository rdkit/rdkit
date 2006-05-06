#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
import time,random

def get12DigitString():
  """ returns a basic 12-digit GUID as a string
  """
  # 9 digits of time:
  t = '%9.0f'%(100*time.time())
  # 3 random digits
  digits = '%03d'%(random.randrange(1000))
  return '%s-%s-%s-%s'%(t[-9:-6],t[-6:-3],t[-3:],digits)


getGUID = get12DigitString
