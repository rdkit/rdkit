#!/usr/bin/env python
# quick script for comparing SVGs in 2 different directories.
import sys
import glob
from pathlib import Path

d1 = '/Users/david/Projects/RDKit/Worktrees/RefactorMolDraw2D/cmake-build-debug/Code/GraphMol/MolDraw2D'
d2 = '/Users/david/Projects/RDKit/gregs_fork/cmake-build-debug/Code/GraphMol/MolDraw2D'

print(f'''<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <meta http-equiv="X-UA-Compatible" content="ie=edge" />
    <title>{d1} vs {d2}</title>
    <link href="index.css" rel="stylesheet" />
  </head>
  <body>
    <table border="1">
      <tr>
        <td></td>
        <td>{d1}</td>
        <td>{d2}</td>
      </tr>''')

inglob = Path(d1) / sys.argv[1]
fns = [Path(fn) for fn in glob.glob(inglob.name)]
fns.sort(key = lambda f: f.stat().st_mtime, reverse=True)

for fp in fns:
    fn = fp.name
    if not fn.endswith('.svg') and not fn.endswith('.png'):
        continue
    fns = fn.replace('.svg', '')
    print(f'''      <tr>
        <td>{fns}</td>
        <td><img src="{d1}/{fn}" alt="{fns}"/></td>
        <td><img src="{d2}/{fn}" alt="{fns}"/></td>
      </tr>''')
    
print('''    </table>
  </body>
</html>''')
