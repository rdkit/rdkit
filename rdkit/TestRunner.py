#
#  Copyright (c) 2003-2010 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import os
import sys
import time

from rdkit import RDConfig

try:
  import subprocess
except ImportError:
  subprocess = None

TEST_FAILED = -1
TEST_PASSED = 0

BUILD_TYPE_ENVVAR = 'RDKIT_BUILD_TYPE'


def isDebugBuild():
  try:
    return os.environ[BUILD_TYPE_ENVVAR] == 'DEBUG'
  except KeyError:
    return False


def RunTest(exeName, args, extras):
  if exeName == "python":
    exeName = RDConfig.pythonExe
  args = args.split(' ')

  startDir = os.getcwd()
  if 'dir' in extras:
    os.chdir(extras['dir'])
  expectedReturn = extras.get('returns', 0)
  if not subprocess:
    raise NotImplementedError('cannot run tests if the subprocess module is not available.')
  else:
    try:
      retVal = subprocess.call([exeName] + list(args))
    except OSError:
      print("Could not find executable: %s." % exeName, file=sys.stderr)
      return TEST_FAILED
  if 'dir' in extras:
    os.chdir(startDir)
  if retVal != expectedReturn:
    return TEST_FAILED
  else:
    return TEST_PASSED


def RunScript(script, doLongTests, verbose):
  # support python 2.7 style -f argument for failfast
  if sys.argv[-1] == '-f':
    # setting environment allows this setting to recursively pass to all child
    # processes
    os.environ['PYTHON_TEST_FAILFAST'] = '1'
  if len(sys.argv) >= 3 and sys.argv[1] == '--testDir':
    os.chdir(sys.argv[2])
  if len(sys.argv) >= 5 and sys.argv[3] == '--buildType':
    os.environ[BUILD_TYPE_ENVVAR] = sys.argv[4].upper()
  # -------------------------------------------------------
  # this is pretty funny.  Whatever directory we started python in
  # will be in the search path, and if we've changed to another
  # directory, that'll be there too.  HOWEVER... the starting
  # directory will be searched first (at least in python2.2), so we
  # need to make sure that '.' is at the front of the search path
  if sys.path[0] != '.':
    sys.path = ['.'] + sys.path
  script = script.split('.py')[0]
  mod = __import__(script)
  try:
    tests = mod.tests
  except AttributeError:
    return [], 0

  longTests = []
  if doLongTests:
    try:
      longTests = mod.longTests
    except AttributeError:
      pass

  failed = []
  for i, entry in enumerate(tests):
    try:
      exeName, args, extras = entry
    except ValueError:
      print('bad entry:', entry)
      sys.exit(-1)
    try:
      res = RunTest(exeName, args, extras)
    except Exception:
      import traceback
      traceback.print_exc()
      res = TEST_FAILED
    if res != TEST_PASSED:
      failed.append((exeName, args, extras))
      # check failfast setting
      if os.environ.get('PYTHON_TEST_FAILFAST', '') == '1':
        # return immediately
        sys.stderr.write("Exiting from %s\n" % str([exeName] + list(args)))
        return failed, i + 1
  for i, (exeName, args, extras) in enumerate(longTests):
    res = RunTest(exeName, args, extras)
    if res != TEST_PASSED:
      failed.append((exeName, args, extras))
      if os.environ.get('PYTHON_TEST_FAILFAST', '') == '1':
        # return immediately
        sys.stderr.write("Exitng from %s\n" % str([exeName] + list(args)))
        return failed, len(tests) + i + 1

  nTests = len(tests) + len(longTests)
  del sys.modules[script]
  if verbose and failed:
    for exeName, args, extras in failed:
      print("!!! TEST FAILURE: ", exeName, args, extras, file=sys.stderr)
  return failed, nTests


def ReportResults(script, failedTests, nTests, runTime, verbose, dest):
  if not nTests:
    dest.write('!-!-!-!-!-!-!-!-!-!-!\n')
    dest.write('\tScript: %s.  No tests run!\n' % (script))
  elif not len(failedTests):
    dest.write('-----------------\n')
    dest.write('\tScript: %s.  Passed %d tests in %.2f seconds\n' % (script, nTests, runTime))
  else:
    dest.write('!-!-!-!-!-!-!-!-!-!-!\n')
    dest.write('\tScript: %s.  Failed %d (of %d) tests in %.2f seconds\n' %
               (script, len(failedTests), nTests, runTime))
    if verbose:
      for exeName, args, extras in failedTests:
        dirName = extras.get('dir', '.')
        dirName = os.path.abspath(dirName)
        dest.write('\t\t(%s): %s %s\n' % (dirName, exeName, args))


if __name__ == '__main__':
  import getopt
  args, extras = getopt.getopt(sys.argv[1:], 'lv')
  doLongTests = 0
  verbose = 1
  for arg, val in args:
    if arg == '-l':
      doLongTests = 1
    elif arg == '-v':
      verbose = 0

  pwd = os.getcwd()
  totNumFailed = 0
  totNumRun = 0
  failures = []
  timeAccum = 0.0
  for script in extras:
    try:
      open(script, 'r')
    except IOError:
      sys.stderr.write('ERROR: Test script %s could not be opened.\n' % (script))
    else:
      dirName = os.path.dirname(script)
      scriptBase = os.path.basename(script)
      if dirName:
        os.chdir(dirName)
      try:
        t1 = time.time()
        failed, nTests = RunScript(scriptBase, doLongTests, verbose)
        t2 = time.time()
      except ImportError:
        import traceback
        traceback.print_exc()
        sys.stderr.write('ERROR: Could not import test script %s\n' % (script))
      else:
        runTime = t2 - t1
        ReportResults(script, failed, nTests, runTime, verbose, sys.stderr)
        timeAccum += runTime
      if dirName:
        os.chdir(pwd)
      if len(extras) > 1:
        totNumFailed += len(failed)
        totNumRun += nTests
        if len(failed):
          failures.append(script)

  if totNumRun > 1:
    sys.stderr.write('\n\n-*-*-*-*-*-*- Test Results Summary -*-*-*-*-*-*-\n')
    sys.stderr.write('\t\tTotal run time: %.2f seconds\n' % (timeAccum))
    if totNumFailed:
      sys.stderr.write('!!!!!---  %d Failures in %d tests  ---!!!!!\n' % (totNumFailed, totNumRun))
      sys.stderr.write('\tModules with failures:\n')
      for failure in failures:
        sys.stderr.write('\t\t%s\n' % failure)
    else:
      sys.stderr.write('  All %d tests (in %d modules) passed\n' % (totNumRun, len(extras)))
  sys.exit(totNumFailed)


class _RedirectStream:
  _stream = None

  def __init__(self, new_target):
    self._new_target = new_target
    # We use a list of old targets to make this CM re-entrant
    self._old_targets = []

  def __enter__(self):
    self._old_targets.append(getattr(sys, self._stream))
    setattr(sys, self._stream, self._new_target)
    return self._new_target

  def __exit__(self, exctype, excinst, exctb):
    setattr(sys, self._stream, self._old_targets.pop())


class redirect_stdout(_RedirectStream):
  """Context manager for temporarily redirecting stdout to another file.

        # How to send help() to stderr
        with redirect_stdout(sys.stderr):
            help(dir)

        # How to write help() to a file
        with open('help.txt', 'w') as f:
            with redirect_stdout(f):
                help(pow)
    """
  _stream = "stdout"


class redirect_stderr(_RedirectStream):
  """Context manager for temporarily redirecting stderr to another file."""
  _stream = "stderr"


class OutputRedirectC:
  """Context manager which uses low-level file descriptors to suppress
  output to stdout/stderr, optionally redirecting to the named file(s).

  Suppress all output
  with Silence():
    <code>

  Redirect stdout to file
  with OutputRedirectC(stdout='output.txt', mode='w'):
    <code>

  Redirect stderr to file
  with OutputRedirectC(stderr='output.txt', mode='a'):
    <code>
  http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
  >>>

  """

  def __init__(self, stdout=os.devnull, stderr=os.devnull, mode='wb'):
    self.outfiles = stdout, stderr
    self.combine = (stdout == stderr)
    self.mode = mode
    self.saved_streams = None
    self.fds = None
    self.saved_fds = None
    self.null_fds = None
    self.null_streams = None

  def __enter__(self):
    # save previous stdout/stderr
    self.saved_streams = saved_streams = sys.__stdout__, sys.__stderr__
    self.fds = fds = [s.fileno() for s in saved_streams]
    self.saved_fds = [os.dup(fd) for fd in fds]
    # flush any pending output
    for s in saved_streams:
      s.flush()

    # open surrogate files
    if self.combine:
      null_streams = [open(self.outfiles[0], self.mode, 0)] * 2
      if self.outfiles[0] != os.devnull:
        # disable buffering so output is merged immediately
        sys.stdout, sys.stderr = [os.fdopen(fd, 'wb', 0) for fd in fds]
    else:
      null_streams = [open(f, self.mode, 0) for f in self.outfiles]
    self.null_fds = null_fds = [s.fileno() for s in null_streams]
    self.null_streams = null_streams

    # overwrite file objects and low-level file descriptors
    for null_fd, fd in zip(null_fds, fds):
      os.dup2(null_fd, fd)

  def __exit__(self, *args):
    # flush any pending output
    for s in self.saved_streams:
      s.flush()
    # restore original streams and file descriptors
    for saved_fd, fd in zip(self.saved_fds, self.fds):
      os.dup2(saved_fd, fd)
    sys.stdout, sys.stderr = self.saved_streams
    # clean up
    for s in self.null_streams:
      s.close()
    for fd in self.saved_fds:
      os.close(fd)
    return False
