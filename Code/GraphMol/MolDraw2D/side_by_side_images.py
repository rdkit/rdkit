#!/usr/bin/env python
# quick script for comparing SVGs in 2 different directories.
import argparse
import glob
from pathlib import Path

parser = argparse.ArgumentParser(description='Make an HTML table for comparing images.')
parser.add_argument('--dir1', required=True, help='Name of first directory.')
parser.add_argument('--dir2', required=True, help='Name of second directory.')
parser.add_argument('--file-glob', required=True, help='glob for files in first directory')
parser.add_argument('--outfile', default='side_by_side.html', help='Name of HTML file')

args = parser.parse_args()

d1 = args.dir1
d2 = args.dir2

if not Path(d1).exists():
    print(f'Directory {d1} missing.')
    exit(1)
if not Path(d2).exists():
    print(f'Directory {d2} missing.')
    exit(1)

with open(args.outfile, 'w') as f:
    f.write(f'''<!DOCTYPE html>
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
      </tr>\n''')

    inglob = Path(d1)
    fns = [Path(fn) for fn in inglob.glob(args.file_glob)]
    fns.sort(key = lambda f: f.stat().st_mtime, reverse=True)

    for fp in fns:
        fn = fp.name
        if not fn.endswith('.svg') and not fn.endswith('.png'):
            continue
        fns = fn.replace('.svg', '')
        f.write(f'''      <tr>
        <td>{fns}</td>
        <td><img src="{d1}/{fn}" alt="{fns}"/></td>
        <td><img src="{d2}/{fn}" alt="{fns}"/></td>
        </tr>\n''')
    
    f.write('''    </table>
  </body>
</html>\n''')
