from __future__ import print_function
"""Emscripten test runner.  Emscripten tests are output as java files and
must be run by node.

Furthermore, if we built with optimization, we must chdir to the location
of the .js file to properly find the optimized .js.mem memory
layouts.
"""
import sys, os, commands
import argparse
parser = argparse.ArgumentParser(description="Enumerate a mess of reactions")
parser.add_argument('search', type=str, nargs="*",
                    help="directories to search for tests")
parser.add_argument('--verbose', 
                    default=False, action='store_true', dest='verbose',
                    help="verbose logging")
args = parser.parse_args()

try:
    path = sys.argv[1]
except IndexError:
    print("== Please supply a path to search for emscripten tests", file=sys.stderr)
    sys.exit(1)

for path in args.search:
    if not os.path.exists(path):
        print("== search path %r does not exist"%path, file=sys.stderr)
        sys.exit(1)
        
for path in args.search:
    failed = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if os.path.splitext(f)[-1] == ".js":
                sys.stdout.write("\t%s "%f)
                sys.stdout.flush()
                lastpath = os.path.abspath(os.curdir)
                os.chdir(root)

                status, output = commands.getstatusoutput("node %s"%f)
                if status:
                    print ("[FAILED]")
                    if args.verbose:
                        print (output, file=sys.stderr)
                    failed.append(f)
                else:
                    print ("[OK]")
                os.chdir(lastpath)

print ("The following tests failed:")
for fail in failed:
    print ("\t%s"%fail, "[FAILED]")
                
