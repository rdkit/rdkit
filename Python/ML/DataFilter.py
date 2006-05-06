# $Id: DataFilter.py 5030 2006-03-02 18:46:51Z glandrum $
#
#  Copyright (C) 2001-2006  greg Landrum
#
#   @@ All Rights Reserved  @@
#
""" command line utility for filtering data

"""

import sys,getopt
from ML.Data import DataUtils
from Dbase import DbConnection

def Usage():
  print 'Usage: DataFilter [-f frac -c col -v value] DBName inTableName outFileBase'
  sys.exit(-1)
  
args,extra = getopt.getopt(sys.argv[1:],'f:c:v:')
frac = 0.5
col = -1
filtVal=0
for arg,val in args:
  if arg == '-f':
    frac = float(val)
  elif arg == '-c':
    col = int(val)
  elif arg == '-v':
    filtVal = int(val)
  else:
    Usage()

if len(extra) < 3:
  Usage()

dbName = extra[0]
inTable = extra[1]
outBase = extra[2]

connect = DbConnection.DbConnect(dbName,inTable)
inD = connect.GetData()
newD,rejD = DataUtils.FilterData(inD,filtVal,frac,col)
inD = None

namesAndTypes = connect.GetColumnNamesAndTypes()

outF = open('%s-keep.csv'%(outBase),'w+')
for name,type in namesAndTypes:
  print >> outF,'%s,'%name.strip(),
print >> outF

for pt in newD:
  for d in pt:
    print >> outF, '%s,'%(str(d)),
  print >> outF
outF.close()

outF = open('%s-rej.csv'%(outBase),'w+')
for name,type in namesAndTypes:
  print >> outF,'%s,'%name.strip(),
print >> outF

for pt in rejD:
  for d in pt:
    print >> outF, '%s,'%(str(d)),
  print >> outF
outF.close()
  
