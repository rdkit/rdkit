try: True
except NameError:
    True = 1
    False = 0

import sys

# don't edit these here: do eg.
# import ClientCookie; ClientCookie.DEBUG_STREAM = True
# DEBUG_STREAM is used for errors and warnings if the logging and warnings
# modules respectively aren't present.
DEBUG_STREAM = sys.stderr
# These work like equivalents from logging.  Use logging direct if you
# have 2.3.
from _Debug import getLogger, NOTSET, INFO, DEBUG
USE_BARE_EXCEPT = True

# Import names so that they can be imported directly from the package, like
# this:
#from ClientCookie import <whatever>
from _ClientCookie import VERSION, __doc__, \
     Cookie, \
     CookiePolicy, DefaultCookiePolicy, \
     CookieJar, FileCookieJar, LoadError
from _LWPCookieJar import LWPCookieJar, lwp_cookie_str
from _MozillaCookieJar import MozillaCookieJar
from _MSIECookieJar import MSIECookieJar
from _BSDDBCookieJar import BSDDBCookieJar, CreateBSDDBCookieJar
#from _MSIEDBCookieJar import MSIEDBCookieJar
try:
    from urllib2 import AbstractHTTPHandler
except ImportError:
    pass
else:
    from ClientCookie._urllib2_support import \
         Request, \
         HTTPHandler, build_opener, install_opener, urlopen, \
         HTTPRedirectHandler
    from ClientCookie._urllib2_support import \
         OpenerDirector, BaseProcessor, \
         HTTPRequestUpgradeProcessor, \
         HTTPEquivProcessor, SeekableProcessor, HTTPCookieProcessor, \
         HTTPRefererProcessor, \
         HTTPRefreshProcessor, HTTPErrorProcessor, \
         HTTPResponseDebugProcessor, HTTPRedirectDebugProcessor

    import httplib
    if hasattr(httplib, 'HTTPS'):
        from ClientCookie._urllib2_support import HTTPSHandler
    del AbstractHTTPHandler, httplib
from _Util import http2time
str2time = http2time
del http2time

del sys
