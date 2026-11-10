"""Check Python tests for leaks before interpreter shutdown destroys modules."""
import atexit
import ctypes
import gc
import os
import sys

try:
    _leak_check = ctypes.CDLL(None).__lsan_do_leak_check
    _leak_check.argtypes = ()
    _leak_check.restype = None
except Exception as error:
    # Python ignores sitecustomize exceptions, so fail explicitly if unavailable.
    print(f"Cannot enable Python leak checking: {error}", file=sys.stderr, flush=True)
    os._exit(1)


# atexit runs callbacks in reverse order: collect first, then scan for leaks.
# Separate callbacks let any frames created by garbage collection unwind before
# LSan scans, avoiding reports for the still-running collector's Python frame.
# The one-shot LSan check replaces the automatic process-exit check. Module
# destruction later in Python/Boost shutdown is outside this check; allocations
# lost during test execution remain detectable.
atexit.register(_leak_check)
atexit.register(gc.collect)
