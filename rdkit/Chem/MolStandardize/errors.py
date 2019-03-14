# -*- coding: utf-8 -*-
"""
molvs.errors
~~~~~~~~~~~~

This module contains exceptions that are raised by MolVS.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


class MolVSError(Exception):
    pass


class StandardizeError(MolVSError):
    pass


class ValidateError(MolVSError):
    pass


class StopValidateError(ValidateError):
    """Called by Validations to stop any further validations from being performed."""
    pass
