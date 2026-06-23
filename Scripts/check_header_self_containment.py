'''

Check that public RDKit C++ headers are self-contained.

A header is "self-contained" if it compiles on its own: a translation unit
whose only contents are `#include <the/header.h>` must succeed. Headers that
silently depend on another header being included first are a frequent source of
build breaks that only surface under a different standard library, compiler, or
include order (for example libstdc++ vs libc++, or when an umbrella header that
used to drag in the missing dependency is removed).

The check compiles each header under Code/ with `-fsyntax-only`, reusing the
compiler and include/define flags from a configured build's
compile_commands.json. Run it against a build tree configured with
-DCMAKE_EXPORT_COMPILE_COMMANDS=ON:

    python Scripts/check_header_self_containment.py --build-dir build

It exits 0 if every checked header is self-contained, 1 otherwise (printing the
first error from each offending header).

Headers that require an optional external dependency (PostgreSQL, Emscripten,
wxWidgets, ...) or that are deliberately not self-contained (fragments meant to
be textually included at a fixed point in a parent header) are skipped; see
SKIP_PATTERNS below.

'''

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from multiprocessing import Pool

# Repo-relative path fragments to skip, each with the reason it is not checkable
# as a standalone translation unit.
SKIP_PATTERNS = [
    # PostgreSQL cartridge: needs the server headers (postgres.h, Datum, ...).
    'Code/PgSQL/',
    # Boost.Python plumbing: needs Python.h / the Python build configuration.
    'Code/RDBoost/PySequenceHolder.h',
    # Emscripten / WebAssembly build only.
    'Code/GraphMol/MolDraw2D/DrawTextJS.h',
    'Code/GraphMol/MolDraw2D/MolDraw2DJS.h',
    'Code/GraphMol/MolDraw2D/DrawTextFTJS.h',
    # wxWidgets renderer.
    'Code/GraphMol/MolDraw2D/MolDraw2Dwx.h',
    # cairo / FreeType renderers (optional external libraries).
    'Code/GraphMol/MolDraw2D/DrawTextCairo.h',
    'Code/GraphMol/MolDraw2D/MolDraw2DCairo.h',
    'Code/GraphMol/MolDraw2D/DrawTextFT.h',
    'Code/GraphMol/MolDraw2D/DrawTextFTCairo.h',
    'Code/GraphMol/MolDraw2D/DrawTextFTSVG.h',
    'Code/GraphMol/MolDraw2D/Qt/',
    # Deliberately not self-contained: textually included at a fixed point inside
    # a parent header.
    'Code/GraphMol/SubstructLibrary/SubstructLibrarySerialization.h',  # end of SubstructLibrary.h
    'Code/GraphMol/FileParsers/MolSupplier.v1API.h',                   # end of MolSupplier.h
    # Alternate RDValue representation, only compiled when explicitly selected.
    'Code/RDGeneral/RDValue-doublemagic.h',
]

# Directory fragments that never hold public, standalone-checkable headers.
SKIP_DIRS = ('/Wrap/', '/Basement/', '/test', '/Test')


def repo_relative(path, source_root):
  return os.path.relpath(os.path.abspath(path), source_root)


def is_python_binding_header(path):
  """A header that pulls in Python.h / Boost.Python cannot stand alone here."""
  with open(path, encoding='utf-8', errors='ignore') as fh:
    text = fh.read()
  return bool(re.search(r'Python\.h|boost/python|RDBoost/Wrap|<numpy/', text))


def headers(source_root):
  code_dir = os.path.join(source_root, 'Code')
  for root, _, files in os.walk(code_dir):
    norm = root.replace(os.sep, '/') + '/'
    if any(frag in norm for frag in SKIP_DIRS):
      continue
    for name in files:
      if not name.endswith('.h'):
        continue
      # Test fixtures (test*.h / Test*.h) are consumed only by their own test
      # translation unit, not by library consumers, so they are out of scope.
      if name.lower().startswith('test'):
        continue
      path = os.path.join(root, name)
      rel = repo_relative(path, source_root).replace(os.sep, '/')
      if any(pat in rel for pat in SKIP_PATTERNS):
        continue
      if is_python_binding_header(path):
        continue
      yield rel


def build_flags(compile_commands):
  """Union the include dirs and (non-export) defines across all TUs.

  Using the union rather than a single TU's flags gives every header the
  superset of include paths it might need, independent of which library it
  happens to live in. Per-target ``*_EXPORTS`` macros are dropped so headers are
  checked as consumers see them, not as the defining library builds them.

  Returns (compiler, std_flag, flags).
  """
  with open(compile_commands) as fh:
    entries = json.load(fh)
  if not entries:
    raise SystemExit('compile_commands.json is empty')

  includes, defines = [], []
  seen_inc, seen_def = set(), set()
  compiler, std_flag = None, '-std=c++20'
  for entry in entries:
    toks = shlex.split(entry.get('command') or ' '.join(entry['arguments']))
    if compiler is None and toks:
      compiler = toks[0]
    i = 0
    while i < len(toks):
      t = toks[i]
      if t.startswith('-std='):
        std_flag = t
      elif t.startswith('-I'):
        if t not in seen_inc:
          seen_inc.add(t); includes.append(t)
      elif t == '-isystem' and i + 1 < len(toks):
        pair = ('-isystem', toks[i + 1])
        if pair not in seen_inc:
          seen_inc.add(pair); includes += list(pair)
        i += 1
      elif t.startswith('-D') and 'EXPORTS' not in t:
        if t not in seen_def:
          seen_def.add(t); defines.append(t)
      i += 1
  return compiler, std_flag, includes + defines


def check_one(args):
  rel, compiler, std_flag, flags = args
  # The TU includes the header by its path relative to Code/, which is on the
  # include path; this mirrors how the rest of the tree includes it.
  include_path = re.sub(r'^Code/', '', rel)
  tf = tempfile.NamedTemporaryFile(suffix='.cpp', delete=False, mode='w')
  tf.write('#include <%s>\n' % include_path)
  tf.close()
  try:
    proc = subprocess.run(
        [compiler, '-fsyntax-only', std_flag, '-Wno-attributes'] + flags + [tf.name],
        capture_output=True, text=True, timeout=180)
  except subprocess.TimeoutExpired:
    os.unlink(tf.name)
    return (rel, 'timed out')
  os.unlink(tf.name)
  if proc.returncode == 0:
    return (rel, None)
  errors = [ln.strip() for ln in proc.stderr.splitlines() if ' error:' in ln]
  return (rel, errors[0] if errors else (proc.stderr.strip().splitlines() or ['?'])[-1])


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--build-dir', default='build',
                      help='Build directory containing compile_commands.json. Default=%(default)s.')
  parser.add_argument('--source-root', default=None,
                      help='RDKit source root. Default: parent of this script.')
  parser.add_argument('--jobs', type=int, default=os.cpu_count() or 4,
                      help='Parallel compile jobs. Default=%(default)s.')
  parser.add_argument('--compiler', default=None,
                      help='Override the compiler taken from compile_commands.json '
                           '(e.g. to check against a different standard library).')
  args = parser.parse_args(argv)

  source_root = args.source_root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  compile_commands = os.path.join(args.build_dir, 'compile_commands.json')
  if not os.path.isfile(compile_commands):
    raise SystemExit(
        'could not find %s; configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON' % compile_commands)

  compiler, std_flag, flags = build_flags(compile_commands)
  if args.compiler:
    compiler = args.compiler
  hdrs = sorted(headers(source_root))
  print('Checking %d headers with %s %s' % (len(hdrs), compiler, std_flag))

  work = [(h, compiler, std_flag, flags) for h in hdrs]
  with Pool(args.jobs) as pool:
    results = pool.map(check_one, work)

  failures = [(rel, err) for rel, err in results if err]
  if failures:
    print('\n%d header(s) are not self-contained:\n' % len(failures))
    for rel, err in sorted(failures):
      print('  %s' % rel)
      print('      %s' % err)
    return 1
  print('All %d headers are self-contained.' % len(hdrs))
  return 0


if __name__ == '__main__':
  sys.exit(main())
