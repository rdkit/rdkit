from __future__ import print_function
import pyRXP
import sys

parser = pyRXP.Parser()
res = parser.parse(open(sys.argv[1],'r').read())
if res:
  print('SUCCESS')
