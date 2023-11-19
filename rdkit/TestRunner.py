#
#  Copyright (c) 2003-2023 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import os
import sys

BUILD_TYPE_ENVVAR = 'RDKIT_BUILD_TYPE'


def isDebugBuild():
  try:
    return os.environ[BUILD_TYPE_ENVVAR] == 'DEBUG'
  except KeyError:
    return False


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
