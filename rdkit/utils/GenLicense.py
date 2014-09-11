#
# Copyright (C) 2002 greg Landrum and Rational Discovery LLC
#
""" generates license files for our primitive license handler

"""
from __future__ import print_function

from rdkit.utils import Licensing
import sha,base64,time,StringIO
import sys


_salt = Licensing._salt

def DataFromTextDate(dText,mods=None):
  """ expected format: day-month-year

  """
  splitD= dText.split('-')
  if len(splitD)!=3:
    sys.stderr.write('ERROR: date format is day-month-year\n')
    sys.exit(0)
  dateComponents = map(int,splitD)
  day,month,year = dateComponents
  if month > 12:
    sys.stderr.write('ERROR: date format is day-month-year\n')
    sys.exit(0)
    
  dVal = int(time.mktime((year,month,day,
                          0,0,0,
                          0,0,0)))
  digest = sha.new(base64.decodestring(_salt))
  digest.update(dText)
  
  digest.update(str(dVal))
  if not mods:
    res = """Expiration_Date: %s
Verification: %s  
    """%(dText,digest.hexdigest())
  else:
    digest.update(mods.upper())
    res = """Expiration_Date: %s
Modules: %s
Verification: %s  
    """%(dText,mods,digest.hexdigest())
  return res

  

if __name__ == '__main__':
  import sys
  d = sys.argv[1]
  if len(sys.argv)>2:
    mods = ','.join([x.strip() for x in sys.argv[2:]])
  else:
    mods = None
  print(DataFromTextDate(d,mods))


  

