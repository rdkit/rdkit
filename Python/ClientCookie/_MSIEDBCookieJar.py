"""Persistent CookieJar based on MS Internet Explorer cookie database.

Copyright 2003 John J Lee <jjl@pobox.com>

This code is free software; you can redistribute it and/or modify it under
the terms of the BSD License (see the file COPYING included with the
distribution).

**********************************************************************
THIS DOESN'T WORK!

It's just a sketch, to check the base class is OK.

**********************************************************************

"""

from ClientCookie import MSIEBase
from _Util import time2netscape

def set_cookie_hdr_from_cookie(cookie):
    params = []
    if cookie.name is not None:
        params.append("%s=%s" % cookie.name, cookie.value)
    else:
        params.append(cookie.name)
    if cookie.expires:
        params.append("expires=" % time2netscape(cookie.expires))
    if cookie.domain_specified:
        params.append("Domain=%s" % cookie.domain)
    if cookie.path_specified:
        params.append("path=%s" % cookie.path)
    if cookie.port_specified:
        if cookie.port is None:
            params.append("Port")
        else:
            params.append("Port=%s" % cookie.port)
    if cookie.secure:
        params.append("secure")
##     if cookie.comment:
##         params.append("Comment=%s" % cookie.comment)
##     if cookie.comment_url:
##         params.append("CommentURL=%s" % cookie.comment_url)
    return "; ".join(params)

class MSIEDBCookieJar(MSIEBase):
    """A CookieJar that relies on MS Internet Explorer's cookie database.

    XXX Require ctypes or write C extension?  win32all probably requires
    latter.

    **********************************************************************
    THIS DOESN'T WORK!

    It's just a sketch, to check the base class is OK.

    **********************************************************************

    MSIEDBCookieJar, unlike MSIECookieJar, keeps no state for itself, but
    relies on the MS Internet Explorer's cookie database.  It uses the win32
    API functions InternetGetCookie() and InternetSetCookie(), from the wininet
    library.

    Note that MSIE itself may impose additional conditions on cookie processing
    on top of that done by CookiePolicy.  For cookie setting, the class tries
    to foil that by providing the request details and Set-Cookie header it
    thinks MSIE wants to see.  For returning cookies to the server, it's up to
    MSIE.

    Note that session cookies ARE NOT written to disk and won't be accessible
    from other processes.  .clear_session_cookies() has no effect.

    .clear_expired_cookies() has no effect: MSIE is responsible for this.

    .clear() will raise NotImplementedError unless all three arguments are
    given.

    """
    def clear_session_cookies(self): pass
    def clear_expired_cookies(self): pass
    def clear(self, domain=None, path=None, name=None):
        if None in [domain, path, name]:
            raise NotImplementedError()
        # XXXX
        url = self._fake_url(domain, path)
        hdr = "%s=; domain=%s; path=%s; max-age=0" % (name, domain, path)
        r = windll.InternetSetCookie(url, None, hdr)
        # XXX return value of InternetSetCookie?
    def _fake_url(self, domain, path):
        # to convince MSIE that Set-Cookie is OK
        return "http://%s%s" % (domain, path)
    def set_cookie(self, cookie):
        # XXXX
        url = self._fake_url(cookie.domain, cookie.path)
        r = windll.InternetSetCookie(
            url, None, set_cookie_hdr_from_cookie(cookie))
        # XXX return value of InternetSetCookie?
    def add_cookie_header(self, request, unverifiable=False):
        # XXXX
        cookie_header = windll.InternetGetCookie(request.get_full_url())
        # XXX return value of InternetGetCookie?
        request.add_unredirected_header(cookie_header)
    def __iter__(self):
        self._load_index_dat()
        return CookieJar.__iter__(self)
    def _cookies_for_domain(self, domain, request, unverifiable):
        #raise NotImplementedError()  # XXXX
        debug("Checking %s for cookies to return", domain)
        if not self.policy.domain_return_ok(domain, request, unverifiable):
            return []

        # XXXX separate out actual loading of cookie data, so only index.dat is
        #  read in ._load_index_dat(), and ._really_load() calls that, then
        #  ._delayload_domain for all domains if not self.delayload.
        #  We then just call ._load_index_dat()
        self._delayload = False
        self._really_load()

        cookies_by_path = self._cookies.get(domain)
        if cookies_by_path is None:
            return []

        cookies = []
        for path in cookies_by_path.keys():
            if not self.policy.path_return_ok(path, request, unverifiable):
                continue
            for name, cookie in cookies_by_path[path].items():
                if not self.policy.return_ok(cookie, request, unverifiable):
                    debug("   not returning cookie")
                    continue
                debug("   it's a match")
                cookies.append(cookie)

        return cookies

