#
# copyright (c) 2003 Rational Discovery LLC
#  All Rights Reserved
#
""" A simple, primitive, shelve-based data storage system for persistence
  on web pages.

"""
import os,shelve,time

def EncodeData(inD,level=5):
  import zlib,base64,urllib
  return urllib.quote_plus(base64.encodestring(zlib.compress(inD,level)))
def DecodeData(inD):
  import zlib,base64,urllib
  return zlib.decompress(base64.decodestring(urllib.unquote_plus(inD)))

_store = {}
def InitStore(ourCookie,fName):
  global _store
  _store[ourCookie] = shelve.open(fName)
  _store[ourCookie]['_lastAccessed']=time.time()
  
def StoreData(ourCookie,what,key):
  global _store
  db = _store[ourCookie]
  db[key] = what
  _store[ourCookie]['_lastModified']=time.time()

def RetrieveData(ourCookie,key):
  global _store
  try:
    db = _store[ourCookie]
    val = db[key]
  except KeyError:
    val = None
  return val
 
def CloseStore(ourCookie):
  try:
    db = _store[ourCookie]
  except:
    pass
  else:
    db.close()
    del _store[ourCookie]

def DeleteStore(ourCookie,fName):

  try:
    os.unlink(fName)
  except OSError:
    pass
