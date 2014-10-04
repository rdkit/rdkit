# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" primitive license handler

  License file format:
  (lines beginning with # are ignored)
    Expiration_Date: <expiration date of license>
    Verification: <verification code>

  Date format: day-month-year
  The verification code is used to ensure that the date
  has not been tampered with
"""
from __future__ import print_function

import sha,base64,time


EXPIRED=-1
BADMODULE=0

class LicenseError(Exception):
  pass

# this is a base64 encoding of the string "RD License Manager"
# it's done this way to provide minimal security (by preventing
#  a simple run of "strings")
_salt='UkQgTGljZW5zZSBNYW5hZ2Vy\n'
#
# Verification strings are constructed using the SHA
#  digest of the message formed from:
#   1) the results of base64 decoding _salt (defined above)
#   2) the text of the Expiration Date
#   3) the string representation of the int form of the
#      time.mktime() date corresponding to the Expiration Date

def _Verify(lines):
  verifier = sha.new(base64.decodestring(_salt))
  inL = lines[0]
  if inL == '':
    raise LicenseError('EOF hit parsing license file')

  if inL.find('Expiration_Date:')==-1:
    raise LicenseError('bad license file format')
  dText = inL.split(':')[-1].strip()
  verifier.update(dText)
  try:
    dateComponents = map(int,dText.split('-'))
  except:
    dateComponents = []
  if len(dateComponents) != 3:
    raise LicenseError('bad date format in license file')
  day,month,year = dateComponents

  pos = 1
  for line in lines:
    if line.find('Modules:')!=-1:
      break
  if line.find('Modules:') != -1:
    modules = line.split(':')[-1].strip().upper().split(',')
    #modules = ','.join([x.strip() for x in modules.split(',')])
  else:
    modules = ''
    
  pos = 1
  inL = lines[pos]
  while pos < len(lines) and inL.find('Verification:')==-1:
    pos += 1
    inL = lines[pos]
  if inL == '':
    raise LicenseError('EOF hit parsing license file')
  if inL.find('Verification:')==-1:
    raise LicenseError('bad license file format')
  vText = inL.split(':')[-1].strip()

  expDate = int(time.mktime((year,month,day,
                             0,0,0,
                             0,0,0)))

  verifier.update(str(expDate))
  verifier.update(','.join(modules))
  if verifier.hexdigest() != vText:
    raise LicenseError('verification of license file failed')

  # ok, the license file has not been tampered with... proceed
  return expDate,modules

def _CheckDate(expDate):
  year,month,day,h,m,s,w,d,dst = time.gmtime()
  newD = int(time.mktime((year,month,day,
                          0,0,0,
                          0,0,0)))
  if expDate > newD:
    return 1
  else:
    return 0
             
  
def CheckLicenseFile(filename):
  try:
    inF = open(filename,'r')
  except IOError:
    raise LicenseError('License file %s could not be opened'%(filename))

  lines = []
  for line in inF.readlines():
    if len(line) and line[0] != '#':
      lines.append(line.strip())
  expDate,modules = _Verify(lines)
  if not _CheckDate(expDate):
    return EXPIRED

def CheckLicenseString(text,checkModule=None):
  lines = text.split('\n')
  expDate,modules = _Verify(lines)
  if not _CheckDate(expDate):
    return EXPIRED
  if checkModule:
    if checkModule.upper() not in modules:
      return BADMODULE
      
  return 1  
if __name__ == '__main__':
  import sys
  if len(sys.argv)>1:
    for nm in sys.argv[1:]:
      print(nm,CheckLicenseFile(nm),CheckLicenseString(open(nm,'r').read()))


  



