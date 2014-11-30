#
#  Copyright (C) 2000  greg Landrum
#
""" contains class REFile, for working with files containing comments

"""
from __future__ import print_function

import re

class REFile:
  """behaves more or less like a normal file, but removes comments

   Any read from the file skips lines beginning with a comment character
     and cleaves off any portion of the line following a comment character

  """
  def close(self):
    """ closes the file

    """
    self.inFile.close()
    
  def rewind(self):
    """ rewinds the file (seeks to the beginning)

    """
    self.inFile.seek(0)
    
  def readlines(self):
    """ reads in all the lines of the file

    """
    res = []
    eofHit = 0
    while not eofHit:
      l = self.inFile.readline()
      if l == '':
        eofHit = 1
      else:
        l = self.expr.split(l)[0]
        eofHit = 0
      while len(l) == 0 and not eofHit:
        l = self.inFile.readline()
        if l == '':
          eofHit = 1
        else:
          l = self.expr.split(l)[0]
      if not eofHit:
        res.append(l)
    return res
  
  def readline(self):
    """ reads in a single line from the file

    """
    l = self.inFile.readline()
    if l == '':
      eofHit = 1
    else:
      l = self.expr.split(l)[0]
      eofHit = 0
    while len(l) == 0 and not eofHit:
      l = self.inFile.readline()
      if l == '':
        eofHit = 1
      else:
        l = self.expr.split(l)[0]
    return l
    
  def __init__(self,fileName,mode='r',commentChar='#|\n'):
    """ Constructor

      **Arguments**

        - fileName: the filename from which to read

        - mode: the mode in which to open the file

        - commentChar: a regexp defining the comment character

    """    
    self.expr = re.compile(commentChar)
    self.fileName = fileName
    self.mode = mode
    self.inFile = open(fileName,mode)

    

if __name__ == '__main__':
  fName = 'retest.txt'
  ref = REFile(fName)
  lines = ref.readlines()
  print('readlines:')
  for i in xrange(len(lines)):
    print('\t%d: %s'%(i,lines[i]))
  ref.rewind()
  print('readline:')
  inStr = ref.readline()
  nRead = 0
  while inStr != '':
    print('\t%d: %s'%(nRead,inStr))
    nRead=nRead+1
    inStr = ref.readline()
  
    
  
