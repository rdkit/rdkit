from cPickle import load as _load
from cPickle import loads as _loads
from cPickle import *

def load(f, **kwargs): return _load(f)
def loads(s, **kwargs): return _loads(s)
