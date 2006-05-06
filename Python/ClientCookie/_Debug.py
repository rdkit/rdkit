import ClientCookie

try:
    import warnings
except ImportError:
    def warn(text):
        ClientCookie.DEBUG_STREAM.write("WARNING: "+text)
else:
    def warn(text):
        warnings.warn(text, stacklevel=2)

try:
    import logging
except:
    NOTSET = None
    INFO = 20
    DEBUG = 10
    class Logger:
        def __init__(self):
            self.level = NOTSET
        def log(self, level, text, *args):
            if args:
                text = text % args
            if self.level is not None and level < self.level:
                ClientCookie.DEBUG_STREAM.write(text+"\n")
        def debug(self, text, *args):
            apply(self.log, (DEBUG, text)+args)
        def info(self, text, *args):
            apply(self.log, (INFO, text)+args)
        def setLevel(self, lvl):
            self.level = lvl
    LOGGER = Logger()
    def getLogger(name): return LOGGER
else:
    from logging import getLogger, INFO, DEBUG, NOTSET
    getLogger("ClientCookie").addHandler(
        logging.StreamHandler(ClientCookie.DEBUG_STREAM))

## def _debug(text, *args):
##     if args:
##         text = text % args
##     ClientCookie.DEBUG_STREAM.write(text+"\n")
