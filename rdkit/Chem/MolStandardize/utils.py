# -*- coding: utf-8 -*-
"""
molvs.utils
~~~~~~~~~~~

This module contains miscellaneous utility functions.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import functools
from itertools import tee

import six


def memoized_property(fget):
    """Decorator to create memoized properties."""
    attr_name = '_{}'.format(fget.__name__)

    @functools.wraps(fget)
    def fget_memoized(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fget(self))
        return getattr(self, attr_name)
    return property(fget_memoized)


def pairwise(iterable):
    """Utility function to iterate in a pairwise fashion."""
    a, b = tee(iterable)
    next(b, None)
    return six.moves.zip(a, b)
