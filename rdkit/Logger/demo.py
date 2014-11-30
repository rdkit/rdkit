# $Id$
#
#       Copyright (c) 2001-2006, Greg Landrum and Rational Discovery LLC,
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

# Demo code for the Logger class
from __future__ import print_function
import Logger

class Foo:
  """ a simple class
  """
  def __init__(self,aVal):
    self.a = aVal
    
  def method1(self,a,b,c='foo'):
    print('method1',a,b,c)
    return 'method1'
  def method2(self):
    print('method2')
    return 'method2'

def demo1():
  l = Logger.Logger(Foo,7)
  l.method1(1,2)
  l.method1(1,2,c='grm')
  l.method2()
  l.method1(7,6,'pizza')
  l.b = 3
  l.method2()
  
  print(l._LoggerGetLog())

  f = Foo(6)
  r = Logger.replay(l._LoggerGetLog(),f)
  print('playback results:',r)
  # f is now in more or less the same state as l... the only differences
  #  will arise because of different arguments passed to the underlying
  #  classes __init__ method
  print(f.b)
  
def demo2():
  # create a Logger which will flush itself:
  l = Logger.Logger(Foo,7,loggerFlushCommand='method2')
  l.method1(23,42)
  l.method1(1,2,'foo')
  print(l._LoggerGetLog())
  # this will blow out the log
  l.method2()
  print(l._LoggerGetLog())

if __name__ == '__main__':
  demo1()
  demo2()
