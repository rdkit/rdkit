"""Integration with Python standard library module urllib2.

Also includes a redirection bugfix, support for parsing HTML HEAD blocks for
the META HTTP-EQUIV tag contents, and following Refresh header redirects.

Copyright 2002-2003 John J Lee <jjl@pobox.com>

This code is free software; you can redistribute it and/or modify it under
the terms of the BSD License (see the file COPYING included with the
distribution).

"""

import copy, time

import ClientCookie
from _ClientCookie import CookieJar, request_host
from _Util import isstringlike, startswith, getheaders
from _Debug import getLogger
info = getLogger("ClientCookie").info

try: True
except NameError:
    True = 1
    False = 0


CHUNK = 1024  # size of chunks fed to HTML HEAD parser, in bytes

def methnames(obj):
    """Return method names of class instance.

    dir(obj) doesn't work across Python versions, this does.

    """
    return methnames_of_instance_as_dict(obj).keys()

def methnames_of_instance_as_dict(inst):
    names = {}
    names.update(methnames_of_class_as_dict(inst.__class__))
    for methname in dir(inst):
        candidate = getattr(inst, methname)
        if callable(candidate):
            names[methname] = None
    return names

def methnames_of_class_as_dict(klass):
    names = {}
    for methname in dir(klass):
        candidate = getattr(klass, methname)
        if callable(candidate):
            names[methname] = None
    for baseclass in klass.__bases__:
        names.update(methnames_of_class_as_dict(baseclass))
    return names

try:
    from urllib2 import AbstractHTTPHandler
except ImportError:
    pass
else:
    import urlparse, urllib2, urllib, httplib#, robotparser
    from urllib2 import URLError, HTTPError
    import types, string, socket, bisect
    from cStringIO import StringIO
    from _Util import response_seek_wrapper
    try:
        import threading
        _threading = threading; del threading
    except ImportError:
        import dummy_threading
        _threading = dummy_threading; del dummy_threading

    # This fixes a bug in urllib2 as of Python 2.1.3 and 2.2.2
    #  (http://www.python.org/sf/549151)
    # 2.2.3 is broken here (my fault!), 2.3 is fixed.
    class HTTPRedirectHandler(urllib2.BaseHandler):
        # maximum number of redirections to any single URL
        # this is needed because of the state that cookies introduce
        max_repeats = 4
        # maximum total number of redirections (regardless of URL) before
        # assuming we're in a loop
        max_redirections = 10

        # Implementation notes:

        # To avoid the server sending us into an infinite loop, the request
        # object needs to track what URLs we have already seen.  Do this by
        # adding a handler-specific attribute to the Request object.  The value
        # of the dict is used to count the number of times the same url has
        # been visited.  This is needed because this isn't necessarily a loop:
        # there is more than one way to redirect (Refresh, 302, 303, 307).

        # Another handler-specific Request attribute, original_url, is used to
        # remember the URL of the original request so that it is possible to
        # decide whether or not RFC 2965 cookies should be turned on during
        # redirect.

        # Always unhandled redirection codes:
        # 300 Multiple Choices: should not handle this here.
        # 304 Not Modified: no need to handle here: only of interest to caches
        #     that do conditional GETs
        # 305 Use Proxy: probably not worth dealing with here
        # 306 Unused: what was this for in the previous versions of protocol??

        def redirect_request(self, newurl, req, fp, code, msg, headers):
            """Return a Request or None in response to a redirect.

            This is called by the http_error_30x methods when a redirection
            response is received.  If a redirection should take place, return a
            new Request to allow http_error_30x to perform the redirect;
            otherwise, return None to indicate that an HTTPError should be
            raised.

            """
            if code in (301, 302, 303, "refresh") or \
                   (code == 307 and not req.has_data()):
                # Strictly (according to RFC 2616), 301 or 302 in response to
                # a POST MUST NOT cause a redirection without confirmation
                # from the user (of urllib2, in this case).  In practice,
                # essentially all clients do redirect in this case, so we do
                # the same.
                return Request(newurl, headers=req.headers)
            else:
                raise HTTPError(req.get_full_url(), code, msg, headers, fp)

        def http_error_302(self, req, fp, code, msg, headers):
            # Some servers (incorrectly) return multiple Location headers
            # (so probably same goes for URI).  Use first header.
            if headers.has_key('location'):
                newurl = getheaders(headers, 'location')[0]
            elif headers.has_key('uri'):
                newurl = getheaders(headers, 'uri')[0]
            else:
                return
            newurl = urlparse.urljoin(req.get_full_url(), newurl)

            # XXX Probably want to forget about the state of the current
            # request, although that might interact poorly with other
            # handlers that also use handler-specific request attributes
            new = self.redirect_request(newurl, req, fp, code, msg, headers)
            if new is None:
                return

            # remember where we started from
            try: new.origin_req_host = req.origin_req_host
            except AttributeError: pass
            new.unverifiable = True

            # loop detection
            # .redirect_dict has a key url if url was previously visited.
            if hasattr(req, 'redirect_dict'):
                visited = new.redirect_dict = req.redirect_dict
                if (visited.get(newurl, 0) >= self.max_repeats or
                    len(visited) >= self.max_redirections):
                    raise HTTPError(req.get_full_url(), code,
                                    self.inf_msg + msg, headers, fp)
            else:
                visited = new.redirect_dict = req.redirect_dict = {}
            visited[newurl] = visited.get(newurl, 0) + 1

            # Don't close the fp until we are sure that we won't use it
            # with HTTPError.  
            fp.read()
            fp.close()

            return self.parent.open(new)

        http_error_301 = http_error_303 = http_error_307 = http_error_302
        http_error_refresh = http_error_302

        inf_msg = "The HTTP server returned a redirect error that would " \
                  "lead to an infinite loop.\n" \
                  "The last 30x error message was:\n"


    class Request(urllib2.Request):
        def __init__(self, url, data=None, headers={}):
            urllib2.Request.__init__(self, url, data, headers)
            self.unredirected_hdrs = {}
            # flag indicating whether the transaction is unverifiable, as
            # defined by RFC 2965
            self.unverifiable = False
            self.origin_req_host = request_host(self)

        def add_unredirected_header(self, key, val):
            """Add a header that will not be added to a redirected request."""
            self.unredirected_hdrs[string.capitalize(key)] = val

        def has_header(self, header_name):
            """True iff request has named header (regular or unredirected)."""
            if (self.headers.has_key(header_name) or
                self.unredirected_hdrs.has_key(header_name)):
                return True
            return False

        def get_header(self, header_name, default=None):
            return self.headers.get(
                header_name,
                self.unredirected_hdrs.get(header_name, default))

        def header_items(self):
            hdrs = self.unredirected_hdrs.copy()
            hdrs.update(self.headers)
            return hdrs.items()

        def __str__(self):
            return "<Request for %s>" % self.get_full_url()

    def print_stack():
        try:
            raise SyntaxError
        except:
            import sys
        import traceback
        traceback.print_stack(sys.exc_traceback.tb_frame.f_back.f_back)

    class BaseProcessor:
        processor_order = 500

        def add_parent(self, parent):
            self.parent = parent
        def close(self):
            self.parent = None
        def __cmp__(self, other):
            if not hasattr(other, "processor_order"):
                return 0
            return cmp(self.processor_order, other.processor_order)
##         def __lt__(self, other):
##             if not hasattr(other, "processor_order"):
##                 return True
##             return self.processor_order < other.processor_order

    class HTTPRequestUpgradeProcessor(BaseProcessor):
        # upgrade Request to class with support for headers that don't get
        # redirected
        processor_order = 0  # before anything else

        def http_request(self, request):
            if not hasattr(request, "add_unredirected_header"):
                newrequest = Request(request._Request__original, request.data,
                                     request.headers)
                # yuck
                try: newrequest.origin_req_host = request.origin_req_host
                except AttributeError: pass
                try: newrequest.unverifiable = request.unverifiable
                except AttributeError: pass
                request = newrequest
            return request

        https_request = http_request

    class HTTPEquivProcessor(BaseProcessor):
        """Append META HTTP-EQUIV headers to regular HTTP headers."""
        def http_response(self, request, response):
            if not hasattr(response, "seek"):
                response = response_seek_wrapper(response)
            # grab HTTP-EQUIV headers and add them to the true HTTP headers
            headers = response.info()
            for hdr, val in parse_head(response):
                # rfc822.Message interprets this as appending, not clobbering
                headers[hdr] = val
            response.seek(0)
            return response

        https_response = http_response

    # XXX ATM this only takes notice of http responses -- probably
    #   should be independent of protocol scheme (http, ftp, etc.)
    class SeekableProcessor(BaseProcessor):
        """Make responses seekable."""

        def http_response(self, request, response):
            if not hasattr(response, "seek"):
                return response_seek_wrapper(response)
            return response

        https_response = http_response

    class HTTPCookieProcessor(BaseProcessor):
        """Handle HTTP cookies.

        Public attributes:

        cookiejar: CookieJar instance

        """
        def __init__(self, cookiejar=None):
            if cookiejar is None:
                cookiejar = CookieJar()
            self.cookiejar = cookiejar

        def http_request(self, request):
            self.cookiejar.add_cookie_header(request)
            return request

        def http_response(self, request, response):
            self.cookiejar.extract_cookies(response, request)
            return response

        https_request = http_request
        https_response = http_response

##     class HTTPRobotRulesProcessor(BaseProcessor):
##         def __init__(self, rfp):
##             self.rfp = rfp
##         def http_request(self, request):
##             host = request.get_host()
##             robots_url = 1
##             return request
##         https_request = http_request

    class HTTPRefererProcessor(BaseProcessor):
        """Add Referer header to requests.

        This only makes sense if you use each RefererProcessor for a single
        chain of requests only (so, for example, if you use a single
        HTTPRefererProcessor to fetch a series of URLs extracted from a single
        page, this will break).

        """
        def __init__(self):
            self.referer = None

        def http_request(self, request):
            if ((self.referer is not None) and
                not request.has_header("Referer")):
                request.add_unredirected_header("Referer", self.referer)
            return request

        def http_response(self, request, response):
            self.referer = response.geturl()
            return response

        https_request = http_request
        https_response = http_response

    class HTTPResponseDebugProcessor(BaseProcessor):
        processor_order = 900  # before redirections, after everything else

        def http_response(self, request, response):
            if not hasattr(response, "seek"):
                response = response_seek_wrapper(response)
            info(response.read())
            info("*****************************************************")
            response.seek(0)
            return response

        https_response = http_response

    class HTTPRedirectDebugProcessor(BaseProcessor):
        def http_request(self, request):
            if hasattr(request, "redirect_dict"):
                info("redirecting to %s", request.get_full_url())
            return request

    class HTTPRefreshProcessor(BaseProcessor):
        """Perform HTTP Refresh redirections.

        Note that if a non-200 HTTP code has occurred (for example, a 30x
        redirect), this processor will do nothing.

        By default, only zero-time Refresh headers are redirected.  Use the
        max_time attribute / constructor argument to allow Refresh with longer
        pauses.  Use the honor_time attribute / constructor argument to control
        whether the requested pause is honoured (with a time.sleep()) or
        skipped in favour of immediate redirection.

        Public attributes:

        max_time: see above
        honor_time: see above

        """
        processor_order = 1000

        def __init__(self, max_time=0, honor_time=True):
            self.max_time = max_time
            self.honor_time = honor_time

        def http_response(self, request, response):
            code, msg, hdrs = response.code, response.msg, response.info()

            if code == 200 and hdrs.has_key("refresh"):
                refresh = getheaders(hdrs, "refresh")[0]
                i = string.find(refresh, ";")
                if i != -1:
                    pause, newurl_spec = refresh[:i], refresh[i+1:]
                    i = string.find(newurl_spec, "=")
                    if i != -1:
                        pause = int(pause)
                        if (self.max_time is None) or (pause <= self.max_time):
                            if pause != 0 and self.honor_time:
                                time.sleep(pause)
                            newurl = newurl_spec[i+1:]
                            # fake a 302 response
                            hdrs["location"] = newurl
                            response = self.parent.error(
                                'http', request, response,
                                "refresh", msg, hdrs)

            return response

        https_response = http_response

    class HTTPErrorProcessor(BaseProcessor):
        """Process HTTP error responses.

        The purpose of this handler is to to allow other response processors a
        look-in by removing the call to parent.error() from
        AbstractHTTPHandler.

        For non-200 error codes, this just passes the job on to the
        Handler.<proto>_error_<code> methods, via the OpenerDirector.error
        method.  Eventually, urllib2.HTTPDefaultErrorHandler will raise an
        HTTPError if no other handler handles the error.

        """
        processor_order = 1000  # after all other processors

        def http_response(self, request, response):
            code, msg, hdrs = response.code, response.msg, response.info()

            if code != 200:
                response = self.parent.error(
                    'http', request, response, code, msg, hdrs)

            return response

        https_response = http_response


    def insort(a, x, lo=0, hi=None, lt=lambda x,y: x<y):
        if hi is None:
            hi = len(a)
        while lo < hi:
            mid = divmod((lo+hi), 2)[0]
            if lt(x, a[mid]): hi = mid
            else: lo = mid+1
        a.insert(lo, x)

    class OpenerDirector(urllib2.OpenerDirector):
        def __init__(self):
            urllib2.OpenerDirector.__init__(self)
            self.process_response = {}
            self.process_request = {}

        def add_handler(self, handler):
            added = False
            for meth in methnames(handler):
                i = string.find(meth, "_")
                protocol = meth[:i]
                condition = meth[i+1:]

                if startswith(condition, "error"):
                    j = string.find(meth[i+1:], "_") + i + 1
                    kind = meth[j+1:]
                    try:
                        kind = int(kind)
                    except ValueError:
                        pass
                    map = self.handle_error.get(protocol, {})
                    self.handle_error[protocol] = map
                    processor = False
                elif (condition == "open" and
                      protocol not in ["do", "proxy"]):  # hack -- see below
                    map = self.handle_open
                    kind = protocol
                    processor = False
                elif (condition in ["response", "request"] and
                      protocol != "redirect"):  # yucky hack
                    # hack above is to fix HTTPRedirectHandler problem, which
                    # appears to above line to be a processor because of the
                    # redirect_request method :-((
                    map = getattr(self, "process_"+condition)
                    kind = protocol
                    processor = True
                else:
                    continue

                if map.has_key(kind):
                    if processor:
                        lt = lambda x,y: x.processor_order < y.processor_order
                    else:
                        lt = lambda x,y: x<y
                    insort(map[kind], handler, lt=lt)
                else:
                    map[kind] = [handler]
                added = True
                continue

            if added:
                # XXX why does self.handlers need to be sorted?
                bisect.insort(self.handlers, handler)
                handler.add_parent(self)

        def _request(self, url_or_req, data):
            if isstringlike(url_or_req):
                req = Request(url_or_req, data)
            else:
                # already a urllib2.Request or ClientCookie.Request instance
                req = url_or_req
                if data is not None:
                    req.add_data(data)
            return req

        def open(self, fullurl, data=None):
            req = self._request(fullurl, data)
            type_ = req.get_type()

            # pre-process request
            # XXX should we allow a Processor to change the type (URL
            #   scheme) of the request?
            meth_name = type_+"_request"
            for processor in self.process_request.get(type_, []):
                meth = getattr(processor, meth_name)
                req = meth(req)

            response = urllib2.OpenerDirector.open(self, req, data)

            # post-process response
            meth_name = type_+"_response"
            for processor in self.process_response.get(type_, []):
                meth = getattr(processor, meth_name)
                response = meth(req, response)

            return response

        def error(self, proto, *args):
            if proto in ['http', 'https']:
                # XXX http[s] protocols are special-cased
                dict = self.handle_error['http'] # https is not different than http
                proto = args[2]  # YUCK!
                meth_name = 'http_error_%s' % proto
                http_err = 1
                orig_args = args
            else:
                dict = self.handle_error
                meth_name = proto + '_error'
                http_err = 0
            args = (dict, proto, meth_name) + args
            result = apply(self._call_chain, args)
            if result:
                return result

            if http_err:
                args = (dict, 'default', 'http_error_default') + orig_args
                return apply(self._call_chain, args)


    # Note the absence of redirect and header-adding code here
    # (AbstractHTTPHandler), and the lack of other clutter that would be
    # here without Processors.
    class AbstractHTTPHandler(urllib2.BaseHandler):
        processor_order = 500

        def __init__(self, debuglevel=0):
            self._debuglevel = debuglevel

        def add_parent(self, parent):
            urllib2.BaseHandler.add_parent(self, parent)

        def set_http_debuglevel(self, level):
            self._debuglevel = level

        def do_request_(self, request):
            host = request.get_host()
            if not host:
                raise URLError('no host given')

            if request.has_data():  # POST
                data = request.get_data()
                if not request.has_header('Content-type'):
                    request.add_unredirected_header(
                        'Content-type',
                        'application/x-www-form-urlencoded')
                if not request.has_header('Content-length'):
                    request.add_unredirected_header(
                        'Content-length', '%d' % len(data))

            scheme, sel = urllib.splittype(request.get_selector())
            sel_host, sel_path = urllib.splithost(sel)
            if not request.has_header('Host'):
                request.add_unredirected_header('Host', sel_host or host)
            for name, value in self.parent.addheaders:
                name = string.capitalize(name)
                if not request.has_header(name):
                    request.add_unredirected_header(name, value)

            return request

        def do_open(self, http_class, req):
            host = req.get_host()
            if not host:
                raise URLError('no host given')

            h = http_class(host) # will parse host:port
            h.set_debuglevel(self._debuglevel)

            if req.has_data():
                h.putrequest('POST', req.get_selector())
            else:
                h.putrequest('GET', req.get_selector())

            for k, v in req.header_items():
                h.putheader(k, v)

            # httplib will attempt to connect() here.  be prepared
            # to convert a socket error to a URLError.
            try:
                h.endheaders()
            except socket.error, err:
                raise URLError(err)
            if req.has_data():
                h.send(req.get_data())

            code, msg, hdrs = h.getreply()
            fp = h.getfile()

            response = urllib.addinfourl(fp, hdrs, req.get_full_url())
            response.code = code
            response.msg = msg

            return response

    # XXX would self.reset() work, instead of raising this exception?
    class EndOfHeadError(Exception): pass
    class AbstractHeadParser:
        # only these elements are allowed in or before HEAD of document
        head_elems = ("html", "head",
                      "title", "base",
                      "script", "style", "meta", "link", "object")

        def __init__(self):
            self.http_equiv = []
        def start_meta(self, attrs):
            http_equiv = content = None
            for key, value in attrs:
                if key == "http-equiv":
                    http_equiv = value
                elif key == "content":
                    content = value
            if http_equiv is not None:
                self.http_equiv.append((http_equiv, content))

        def end_head(self):
            raise EndOfHeadError()

    # use HTMLParser if we have it (it does XHTML), htmllib otherwise
    try:
        import HTMLParser
    except ImportError:
        import htmllib, formatter
        class HeadParser(AbstractHeadParser, htmllib.HTMLParser):
            def __init__(self):
                htmllib.HTMLParser.__init__(self, formatter.NullFormatter())
                AbstractHeadParser.__init__(self)

            def handle_starttag(self, tag, method, attrs):
                if tag in self.head_elems:
                    method(attrs)
                else:
                    raise EndOfHeadError()

            def handle_endtag(self, tag, method):
                if tag in self.head_elems:
                    method()
                else:
                    raise EndOfHeadError()

        HEAD_PARSER_CLASS = HeadParser
    else:
        class XHTMLCompatibleHeadParser(AbstractHeadParser,
                                        HTMLParser.HTMLParser):
            def __init__(self):
                HTMLParser.HTMLParser.__init__(self)
                AbstractHeadParser.__init__(self)

            def handle_starttag(self, tag, attrs):
                if tag not in self.head_elems:
                    raise EndOfHeadError()
                try:
                    method = getattr(self, 'start_' + tag)
                except AttributeError:
                    try:
                        method = getattr(self, 'do_' + tag)
                    except AttributeError:
                        pass # unknown tag
                    else:
                        method(attrs)
                else:
                    method(attrs)

            def handle_endtag(self, tag):
                if tag not in self.head_elems:
                    raise EndOfHeadError()
                try:
                    method = getattr(self, 'end_' + tag)
                except AttributeError:
                    pass # unknown tag
                else:
                    method()

            # handle_charref, handle_entityref and default entitydefs are taken
            # from sgmllib
            def handle_charref(self, name):
                try:
                    n = int(name)
                except ValueError:
                    self.unknown_charref(name)
                    return
                if not 0 <= n <= 255:
                    self.unknown_charref(name)
                    return
                self.handle_data(chr(n))

            # Definition of entities -- derived classes may override
            entitydefs = \
                    {'lt': '<', 'gt': '>', 'amp': '&', 'quot': '"', 'apos': '\''}

            def handle_entityref(self, name):
                table = self.entitydefs
                if name in table:
                    self.handle_data(table[name])
                else:
                    self.unknown_entityref(name)
                    return

            def unknown_entityref(self, ref):
                self.handle_data("&%s;" % ref)

            def unknown_charref(self, ref):
                self.handle_data("&#%s;" % ref)

        HEAD_PARSER_CLASS = XHTMLCompatibleHeadParser

    def parse_head(fileobj):
        """Return a list of key, value pairs."""
        hp = HEAD_PARSER_CLASS()
        while 1:
            data = fileobj.read(CHUNK)
            try:
                hp.feed(data)
            except EndOfHeadError:
                break
            if len(data) != CHUNK:
                # this should only happen if there is no HTML body, or if
                # CHUNK is big
                break
        return hp.http_equiv


    class HTTPHandler(AbstractHTTPHandler):
        def http_open(self, req):
            return self.do_open(httplib.HTTP, req)

        http_request = AbstractHTTPHandler.do_request_

    if hasattr(httplib, 'HTTPS'):
        class HTTPSHandler(AbstractHTTPHandler):
            def https_open(self, req):
                return self.do_open(httplib.HTTPS, req)

            https_request = AbstractHTTPHandler.do_request_

    def build_opener(*handlers):
        """Create an opener object from a list of handlers and processors.

        The opener will use several default handlers and processors, including
        support for HTTP and FTP.

        If any of the handlers passed as arguments are subclasses of the
        default handlers, the default handlers will not be used.

        """
        opener = OpenerDirector()
        default_classes = [
            # handlers
            urllib2.ProxyHandler,
            urllib2.UnknownHandler,
            HTTPHandler,  # from this module (derived from new AbstractHTTPHandler)
            urllib2.HTTPDefaultErrorHandler,
            HTTPRedirectHandler,  # from this module (bugfixed)
            urllib2.FTPHandler,
            urllib2.FileHandler,
            # processors
            HTTPRequestUpgradeProcessor,
            #HTTPEquivProcessor,
            #SeekableProcessor,
            HTTPCookieProcessor,
            #HTTPRefererProcessor,
            #HTTPRefreshProcessor,
            HTTPErrorProcessor
            ]
        if hasattr(httplib, 'HTTPS'):
            default_classes.append(HTTPSHandler)
        skip = []
        for klass in default_classes:
            for check in handlers:
                if type(check) == types.ClassType:
                    if issubclass(check, klass):
                        skip.append(klass)
                elif type(check) == types.InstanceType:
                    if isinstance(check, klass):
                        skip.append(klass)
        for klass in skip:
            default_classes.remove(klass)

        for klass in default_classes:
            opener.add_handler(klass())
        for h in handlers:
            if type(h) == types.ClassType:
                h = h()
            opener.add_handler(h)

        return opener


    _opener = None
    urlopen_lock = _threading.Lock()
    def urlopen(url, data=None):
        global _opener
        if _opener is None:
            urlopen_lock.acquire()
            try:
                if _opener is None:
                    _opener = build_opener()
            finally:
                urlopen_lock.release()
        return _opener.open(url, data)

    def install_opener(opener):
        global _opener
        _opener = opener
