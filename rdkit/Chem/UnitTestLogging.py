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
import sys
import time
import unittest

from rdkit import Chem


class CaptureStream:
  """Helper class that captures output to a file descriptor"""
  def __init__(self, fd):
    self.fd  = fd
    self.old = os.dup(fd)
    self.pr, self.pw = os.pipe()
    os.dup2(self.pw, self.fd)

  def release(self):
    os.dup2(self.old, self.fd)
    os.close(self.pw)

    result = os.read(self.pr, 1024)
    os.close(self.old)
    os.close(self.pr)

    return result.decode('utf-8')


class CapturePython:
  """Helper class that captures Python output"""
  def __init__(self, stream):
    self.stream = stream
    self.new    = io.StringIO()
    self.old    = getattr(sys, self.stream)
    setattr(sys, self.stream, self.new)

  def release(self):
    setattr(sys, self.stream, self.old)
    return self.new.getvalue()


class CaptureLogger(logging.Handler):
  """Helper class that captures Python logger output"""
  def __init__(self, module=None):
    super(CaptureLogger, self).__init__(level=logging.DEBUG)
    self.logs = collections.defaultdict(io.StringIO)

    self.logger = logging.getLogger(module)
    self.oldlevel = self.logger.level
    self.logger.setLevel(logging.DEBUG)
    self.logger.addHandler(self)

  def emit(self, record):
    key = record.levelname
    val = self.format(record)
    self.logs[key].write(val)

  def release(self):
    self.logger.removeHandler(self)
    self.logger.setLevel(self.oldlevel)
    return {k:s.getvalue() for k, s in self.logs.items()}


class CaptureOutput:
  """Helper class that captures all output"""
  def __init__(self):
    self.captured = {}

  def __enter__(self):
    self.osout = CaptureStream(1)
    self.oserr = CaptureStream(2)
    self.pyout = CapturePython('stdout')
    self.pyerr = CapturePython('stderr')
    self.pylog = CaptureLogger('rdkit')
    return self.captured

  def __exit__(self, x, y, z):
    for key, output in self.pylog.release().items():
      self.captured[key] = output

    pyout = self.pyout.release()
    if pyout:
      self.captured['sys.stdout'] = pyout

    pyerr = self.pyerr.release()
    if pyerr:
      self.captured['sys.stderr'] = pyerr

    osout = self.osout.release()
    if osout:
      self.captured['std::cout'] = osout

    oserr = self.oserr.release()
    if oserr:
      self.captured['std::cerr'] = oserr


def timestamp(message):
  ts = time.time()
  if ts % 1 > 0.995:
    # Avoid failures when seconds roll over:
    time.sleep(0.006)
    ts = time.time()

  dt = datetime.datetime.fromtimestamp(ts)
  return dt.strftime('[%H:%M:%S] ') + message

def expect_debug(message):
  expect = timestamp(message)
  Chem.LogDebugMsg(message)
  return expect

def expect_info(message):
  expect = timestamp(message)
  Chem.LogInfoMsg(message)
  return expect

def expect_warning(message):
  expect = timestamp(message)
  Chem.LogWarningMsg(message)
  return expect

def expect_error(message):
  expect = timestamp(message)
  Chem.LogErrorMsg(message)
  return expect


class TestLogToCppStreams(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    Chem.LogToCppStreams()

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


class TestLogToPythonLogger(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    Chem.LogToPythonLogger()

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


class TestLogToPythonStderr(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    Chem.LogToPythonStderr()

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


class TestWrapLogs(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    Chem.LogToCppStreams()
    Chem.WrapLogs()

  def testDebug(self):
    with CaptureOutput() as captured:
      expect = expect_debug('debug') + '\n'
    self.assertEqual(captured, {
      'sys.stderr': 'RDKit DEBUG: ' + expect,
      'std::cerr': expect
    })

  def testInfo(self):
    with CaptureOutput() as captured:
      expect = expect_info('info') + '\n'
    self.assertEqual(captured, {
      'sys.stderr': 'RDKit INFO: ' + expect,
      'std::cout': expect
    })

  def testWarning(self):
    with CaptureOutput() as captured:
      expect = expect_warning('warning') + '\n'
    self.assertEqual(captured, {
      'sys.stderr': 'RDKit WARNING: ' + expect,
      'std::cerr': expect
    })

  def testError(self):
    with CaptureOutput() as captured:
      expect = expect_error('error') + '\n'
    self.assertEqual(captured, {
      'sys.stderr': 'RDKit ERROR: ' + expect,
      'std::cerr': expect
    })


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
