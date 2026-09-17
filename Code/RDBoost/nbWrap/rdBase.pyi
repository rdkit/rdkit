"""Module containing basic definitions for wrapped C++ code"""



rdkitVersion: str = '2026.09.1pre'

boostVersion: str = ''

rdkitBuild: str = 'Darwin|25.6.0|UNIX|AppleClang|64-bit'

def LogToCppStreams() -> None:
    """Initialize RDKit logs with C++ streams"""

def LogToPythonLogger() -> None:
    """Initialize RDKit logs with Python's logging module"""

def LogToPythonStderr() -> None:
    """Initialize RDKit logs with Python's stderr stream"""

def WrapLogs() -> None:
    """Tee the RDKit logs to Python's stderr stream"""

def EnableLog(spec: str) -> None: ...

def DisableLog(spec: str) -> None: ...

def LogStatus() -> str: ...

def LogDebugMsg(msg: str) -> None:
    """Log a message to the RDKit debug logs"""

def LogInfoMsg(msg: str) -> None:
    """Log a message to the RDKit info logs"""

def LogWarningMsg(msg: str) -> None:
    """Log a message to the RDKit warning logs"""

def LogErrorMsg(msg: str) -> None:
    """Log a message to the RDKit error logs"""

def LogMessage(spec: str, msg: str) -> None:
    """Log a message to any rdApp.* log"""

def AttachFileToLog(spec: str, filename: str, delay: int = 100) -> None:
    """Causes the log to write to a file"""

def SeedRandomNumberGenerator(seed: int) -> None:
    """
    Provides a seed to the standard C random number generator
    This does not affect pure Python code, but is relevant to some of the RDKit C++ components.
    """

class BlockLogs:
    """
    Temporarily block logs from outputting while this instance is in scope.
    """

    def __init__(self) -> None: ...

    def __enter__(self) -> BlockLogs: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> None: ...

class CaptureErrorLog:
    """
    Captures messages from rdErrorLog while this instance is in scope.
          Can be used as a context manager. The ``messages`` property is
          accessible both inside the context and after it exits.
          Nesting is supported: inner captures shadow outer ones.

          Example::

            with rdBase.CaptureErrorLog() as capture:
                rdkit_function_that_may_fail()
            print(capture.messages)
    """

    def __init__(self) -> None: ...

    def __enter__(self) -> CaptureErrorLog: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> None: ...

    @property
    def messages(self) -> str:
        """Messages captured from rdErrorLog."""
