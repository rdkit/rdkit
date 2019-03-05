from pickle import load as _load
from pickle import loads as _loads
from pickle import *


def load(f, **kwargs):
  return _load(f)


def loads(s, **kwargs):
  return _loads(s)
