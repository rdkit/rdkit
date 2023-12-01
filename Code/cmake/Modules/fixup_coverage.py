"""This file searches the source directory for BISON and FLEX files that
the coverage tool mistakenly finds in the build tree.
It replaces the paths with the ones from the source tree
n.b. if a file with the same name (i.e. sln.yy) is found twice
 in the source tree, this will break"""

import os
import sys

source_dir, info_file = sys.argv[1:3]
print(source_dir, info_file)

paths = {}
for root, dir, files in os.walk(source_dir):
  for f in files:
    paths[f] = paths.get(f, []) + [os.path.join(root, f)]

lines = open(info_file).readlines()

newlines = []
for line in lines:
  if "SF:" in line:
    fn = line.split("SF:")[-1].strip()
    if not os.path.exists(fn):
      print("Does not exist:", fn.strip())
      head, rest = os.path.split(fn)
      potential = paths[rest]
      if len(potential) == 1:
        line = "SF:" + potential[0]
      else:
        raise NotImplementedError('asdf')
  newlines.append(line)

open(info_file, 'w').write("\n".join(newlines))
