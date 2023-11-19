#
#  Copyright (C) 2022 Kevin Burk
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
unit testing code for various logging schemes
"""

import collections
import datetime
import io
import logging
import os
import re
import sys
import tempfile
import threading
import time
import unittest

import rdkit
from rdkit import rdBase


class CaptureStream:
  """Helper class that captures output to a file descriptor"""

  def __init__(self, fd):
    self.fd = fd
    self.old = os.dup(self.fd)
    self.tmp = tempfile.TemporaryFile()
    os.dup2(self.tmp.fileno(), self.fd)

  def release(self):
    self.tmp.seek(0)
    result = self.tmp.read()

    os.dup2(self.old, self.fd)
    os.close(self.old)
    self.tmp.close()

    return result.decode('utf-8')


class CapturePython:
  """Helper class that captures Python output"""

  def __init__(self, stream):
    self.stream = stream
    self.new = io.StringIO()
    self.old = getattr(sys, self.stream)
    setattr(sys, self.stream, self.new)

  def release(self):
    setattr(sys, self.stream, self.old)
    return self.new.getvalue()


class CaptureLogger(logging.Handler):
  """Helper class that captures Python logger output"""

  def __init__(self, module=None):
    super(CaptureLogger, self).__init__(level=logging.DEBUG)
    self.logs = collections.defaultdict(str)

    self.devnull = open(os.devnull, 'w')
    rdkit.log_handler.setStream(self.devnull)
    rdkit.logger.addHandler(self)

  def handle(self, record):
    key = record.levelname
    val = self.format(record)
    self.logs[key] += val
    return False

  def release(self):
    rdkit.log_handler.setStream(sys.stderr)
    rdkit.logger.removeHandler(self)
    self.devnull.close()
    return self.logs


class CaptureOutput:
  """Helper class that captures all output"""
  timeexpr = re.compile(r'^\[.*?\]')
  timereplacement = '[timestamp]'

  def __init__(self):
    self.captured = {}

  def __enter__(self):
    self.osout = CaptureStream(1)
    self.oserr = CaptureStream(2)
    self.pyout = CapturePython('stdout')
    self.pyerr = CapturePython('stderr')
    self.pylog = CaptureLogger()
    return self.captured

  def __exit__(self, x, y, z):
    for key, output in self.pylog.release().items():
      self.captured[key] = self.timeexpr.sub(self.timereplacement, output)

    pyout = self.pyout.release()
    if pyout:
      self.captured['sys.stdout'] = self.timeexpr.sub(self.timereplacement, pyout)

    pyerr = self.pyerr.release()
    if pyerr:
      self.captured['sys.stderr'] = self.timeexpr.sub(self.timereplacement, pyerr)

    osout = self.osout.release()
    if osout:
      self.captured['std::cout'] = self.timeexpr.sub(self.timereplacement, osout)

    oserr = self.oserr.release()
    if oserr:
      self.captured['std::cerr'] = self.timeexpr.sub(self.timereplacement, oserr)


# Helpers for the non-threaded tests:
def timestamp(message):
  # using actual timestamps is asking for failures during CI
  return f'{CaptureOutput.timereplacement} {message}'


def expect_debug(message):
  expect = timestamp(message)
  rdBase.LogDebugMsg(message)
  return expect


def expect_info(message):
  expect = timestamp(message)
  rdBase.LogInfoMsg(message)
  return expect


def expect_warning(message):
  expect = timestamp(message)
  rdBase.LogWarningMsg(message)
  return expect


def expect_error(message):
  expect = timestamp(message)
  rdBase.LogErrorMsg(message)
  return expect


# Helpers for the threaded tests:
nthreads = 5
nlogs = 50

logger = logging.getLogger("rdkit")
logger.setLevel(logging.DEBUG)


def go(func, *args):
  thread = threading.Thread(target=func, args=args)
  thread.start()
  return thread


def LogDebugs(nlogs, t=1):
  for i in range(1, nlogs + 1):
    rdBase.LogDebugMsg("Debug %d.%d: My dog has fleas!" % (t, i))


def LogInfos(nlogs, t=1):
  for i in range(1, nlogs + 1):
    rdBase.LogInfoMsg("Info %d.%d: Everything is fine." % (t, i))


def LogWarnings(nlogs, t=1):
  for i in range(1, nlogs + 1):
    rdBase.LogWarningMsg("Warning %d.%d: Every good boy does fine." % (t, i))


def LogErrors(nlogs, t=1):
  for i in range(1, nlogs + 1):
    rdBase.LogErrorMsg("Error %d.%d: Intruder detected!" % (t, i))


def LogAllLevels(nlogs, t=1):
  for i in range(1, nlogs + 1):
    rdBase.LogDebugMsg("Debug %d.%d: Headin' out..." % (t, i))
    rdBase.LogInfoMsg("Info %d.%d: There is no cow level." % (t, i))
    rdBase.LogWarningMsg("Warning %d.%d: Nuclear launch detected!" % (t, i))
    rdBase.LogErrorMsg("Error %d.%d: We require more vespene gas." % (t, i))


def RunOneThreadPerLevel(nthreads):
  threads = []
  for i in range(1, nthreads + 1):
    threads.append(go(LogDebugs, nlogs, i))
    threads.append(go(LogInfos, nlogs, i))
    threads.append(go(LogErrors, nlogs, i))
    threads.append(go(LogWarnings, nlogs, i))
  for t in threads:
    t.join()


def RunManyThreadsPerLevel(nthreads):
  threads = []
  for i in range(1, nthreads + 1):
    threads.append(go(LogAllLevels, nlogs, i))
  for t in threads:
    t.join()


class TestLogToCppStreams(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    rdBase.LogToCppStreams()
    rdBase.EnableLog('rdApp.debug')
    rdBase.EnableLog('rdApp.info')

  def testDebug(self):
    with CaptureOutput() as captured:
      expect = expect_debug('debug') + '\n'
    self.assertEqual(captured, {'std::cerr': expect})

  def testInfo(self):
    with CaptureOutput() as captured:
      expect = expect_info('info') + '\n'
    self.assertEqual(captured, {'std::cout': expect})

  def testWarning(self):
    with CaptureOutput() as captured:
      expect = expect_warning('warning') + '\n'
    self.assertEqual(captured, {'std::cerr': expect})

  def testError(self):
    with CaptureOutput() as captured:
      expect = expect_error('error') + '\n'
    self.assertEqual(captured, {'std::cerr': expect})

  def testSynchronous(self):
    with CaptureOutput() as captured:
      LogAllLevels(nlogs)
    cout = captured['std::cout']
    cerr = captured['std::cerr']
    self.assertEqual(cerr.count('Debug'), nlogs)
    self.assertEqual(cout.count('Info'), nlogs)
    self.assertEqual(cerr.count('Warning'), nlogs)
    self.assertEqual(cerr.count('Error'), nlogs)

  def testAsynchronous1(self):
    with CaptureOutput() as captured:
      RunOneThreadPerLevel(nthreads)
    cout = captured['std::cout']
    cerr = captured['std::cerr']
    self.assertEqual(cerr.count('Debug'), nthreads * nlogs)
    self.assertEqual(cout.count('Info'), nthreads * nlogs)
    self.assertEqual(cerr.count('Warning'), nthreads * nlogs)
    self.assertEqual(cerr.count('Error'), nthreads * nlogs)

  def testAsynchronous2(self):
    with CaptureOutput() as captured:
      RunManyThreadsPerLevel(nthreads)
    cout = captured['std::cout']
    cerr = captured['std::cerr']
    self.assertEqual(cerr.count('Debug'), nthreads * nlogs)
    self.assertEqual(cout.count('Info'), nthreads * nlogs)
    self.assertEqual(cerr.count('Warning'), nthreads * nlogs)
    self.assertEqual(cerr.count('Error'), nthreads * nlogs)


class TestLogToPythonLogger(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    rdBase.LogToPythonLogger()

  def testDebug(self):
    with CaptureOutput() as captured:
      expect = expect_debug('debug')
    self.assertEqual(captured, {'DEBUG': expect})

  def testInfo(self):
    with CaptureOutput() as captured:
      expect = expect_info('info')
    self.assertEqual(captured, {'INFO': expect})

  def testWarning(self):
    with CaptureOutput() as captured:
      expect = expect_warning('warning')
    self.assertEqual(captured, {'WARNING': expect})

  def testError(self):
    with CaptureOutput() as captured:
      expect = expect_error('error')
    self.assertEqual(captured, {'ERROR': expect})

  def testSynchronous(self):
    with CaptureOutput() as captured:
      LogAllLevels(nlogs)
    self.assertEqual(captured['DEBUG'].count('Debug'), nlogs)
    self.assertEqual(captured['INFO'].count('Info'), nlogs)
    self.assertEqual(captured['WARNING'].count('Warning'), nlogs)
    self.assertEqual(captured['ERROR'].count('Error'), nlogs)

  def testAsynchronous1(self):
    with CaptureOutput() as captured:
      RunOneThreadPerLevel(nthreads)
    self.assertEqual(captured['DEBUG'].count('Debug'), nthreads * nlogs)
    self.assertEqual(captured['INFO'].count('Info'), nthreads * nlogs)
    self.assertEqual(captured['WARNING'].count('Warning'), nthreads * nlogs)
    self.assertEqual(captured['ERROR'].count('Error'), nthreads * nlogs)

  def testAsynchronous2(self):
    with CaptureOutput() as captured:
      RunManyThreadsPerLevel(nthreads)
    self.assertEqual(captured['DEBUG'].count('Debug'), nthreads * nlogs)
    self.assertEqual(captured['INFO'].count('Info'), nthreads * nlogs)
    self.assertEqual(captured['WARNING'].count('Warning'), nthreads * nlogs)
    self.assertEqual(captured['ERROR'].count('Error'), nthreads * nlogs)


class TestLogToPythonStderr(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    rdBase.LogToPythonStderr()

  def testDebug(self):
    with CaptureOutput() as captured:
      expect = expect_debug('debug') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect})

  def testInfo(self):
    with CaptureOutput() as captured:
      expect = expect_info('info') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect})

  def testWarning(self):
    with CaptureOutput() as captured:
      expect = expect_warning('warning') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect})

  def testError(self):
    with CaptureOutput() as captured:
      expect = expect_error('error') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect})

  def testSynchronous(self):
    with CaptureOutput() as captured:
      LogAllLevels(nlogs)
    output = captured['sys.stderr']
    self.assertEqual(output.count('Debug'), nlogs)
    self.assertEqual(output.count('Info'), nlogs)
    self.assertEqual(output.count('Warning'), nlogs)
    self.assertEqual(output.count('Error'), nlogs)

  def testAsynchronous1(self):
    with CaptureOutput() as captured:
      RunOneThreadPerLevel(nthreads)
    output = captured['sys.stderr']
    self.assertEqual(output.count('Debug'), nthreads * nlogs)
    self.assertEqual(output.count('Info'), nthreads * nlogs)
    self.assertEqual(output.count('Warning'), nthreads * nlogs)
    self.assertEqual(output.count('Error'), nthreads * nlogs)

  def testAsynchronous2(self):
    with CaptureOutput() as captured:
      RunManyThreadsPerLevel(nthreads)
    output = captured['sys.stderr']
    self.assertEqual(output.count('Debug'), nthreads * nlogs)
    self.assertEqual(output.count('Info'), nthreads * nlogs)
    self.assertEqual(output.count('Warning'), nthreads * nlogs)
    self.assertEqual(output.count('Error'), nthreads * nlogs)


class TestWrapLogs(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    rdBase.LogToCppStreams()
    rdBase.WrapLogs()
    rdBase.EnableLog('rdApp.debug')
    rdBase.EnableLog('rdApp.info')

  def testDebug(self):
    with CaptureOutput() as captured:
      expect = expect_debug('debug') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect, 'std::cerr': expect})

  def testInfo(self):
    with CaptureOutput() as captured:
      expect = expect_info('info') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect, 'std::cout': expect})

  def testWarning(self):
    with CaptureOutput() as captured:
      expect = expect_warning('warning') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect, 'std::cerr': expect})

  def testError(self):
    with CaptureOutput() as captured:
      expect = expect_error('error') + '\n'
    self.assertEqual(captured, {'sys.stderr': expect, 'std::cerr': expect})


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
