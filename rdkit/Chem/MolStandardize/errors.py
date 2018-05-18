# -*- coding: utf-8 -*-
"""
molvs.errors
~~~~~~~~~~~~

This module contains exceptions that are raised by MolVS.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division


class MolVSError(Exception):
    pass


class StandardizeError(MolVSError):
    pass


class ValidateError(MolVSError):
    pass


class StopValidateError(ValidateError):
    """Called by Validations to stop any further validations from being performed."""
    pass
