import os
import string

for fname in filter(lambda x: '-' in x, os.listdir(os.curdir)):
  newname = string.replace(fname, '-medium-', '-')
  if newname != fname:
    print(newname)
    os.rename(fname, newname)
