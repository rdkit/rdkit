"""Persistent CookieJar based on bsddb standard library module.

Copyright 2003 John J Lee <jjl@pobox.com>

This code is free software; you can redistribute it and/or modify it under
the terms of the BSD License (see the file COPYING included with the
distribution).

**********************************************************************
THIS IS NOT FULLY TESTED!
**********************************************************************

"""

from _ClientCookie import CookieJar, MappingIterator
from _Debug import getLogger
debug = getLogger("ClientCookie").debug

import bsddb
import cPickle
pickle = cPickle
del cPickle

try: StopIteration
except NameError:
    from _ClientCookie import StopIteration

def CreateBSDDBCookieJar(filename, policy=None):
    """Return a BSDDBCookieJar given a BSDDB filename.

    Use this unless rather than directly using the BSDDBCookieJar constructor
    unless you know what you're doing.

    filename: filename for sleepycat BSDDB database; if the file doesn't exist,
     it will be created; otherwise, it will be opened

    **********************************************************************
    BSDDBCookieJar IS NOT FULLY TESTED!
    **********************************************************************

    """
    db = bsddb.db.DB()
    db.open(filename, bsddb.db.DB_HASH, bsddb.db.DB_CREATE, 0666)
    return BSDDBCookieJar(policy, db)

class BSDDBIterator:
    # XXXX should this use thread lock?
    def __init__(self, cursor):
        iterator = None
        self._c = cursor
        self._i = iterator
    def __iter__(self): return self
    def close(self):
        if self._c is not None:
            self._c.close()
        self._c = self._i = self.next = self.__iter__ = None
    def next(self):
        while 1:
            if self._i is None:
                item = self._c.next()
                if item is None:
                    self.close()
                    raise StopIteration()
                domain, data = item
                self._i = MappingIterator(pickle.loads(data))
            try:
                return self._i.next()
            except StopIteration:
                self._i = None
                continue
    def __del__(self):
        # XXXX will this work?
        self.close()

class BSDDBCookieJar(CookieJar):
    """CookieJar based on a BSDDB database, using the standard bsddb module.

    You should use CreateBSDDBCookieJar instead of the constructor, unless you
    know what you're doing.

    Note that session cookies ARE stored in the database (marked as session
    cookies), and will be written to disk if the database is file-based.  In
    order to clear session cookies at the end of a session, you must call
    .clear_session_cookies().

    Call the .close() method after you've finished using an instance of this
    class.

    **********************************************************************
    THIS IS NOT FULLY TESTED!
    **********************************************************************

    """
    # XXX
    # use transactions to make multiple reader processes possible
    def __init__(self, policy=None, db=None):
        CookieJar.__init__(self, policy)
        del self._cookies
        if db is None:
            db = bsddb.db.DB()
        self._db = db
    def close(self):
        self._db.close()
    def __del__(self):
        # XXXX will this work?
        self.close()
    def clear(self, domain=None, path=None, name=None):
        if name is not None:
            if (domain is None) or (path is None):
                raise ValueError(
                    "domain and path must be given to remove a cookie by name")
        elif path is not None:
            if domain is None:
                raise ValueError(
                    "domain must be given to remove cookies by path")

        db = self._db
        self._cookies_lock.acquire()
        try:
            if domain is not None:
                data = db.get(domain)
                if data is not None:
                    if path is name is None:
                        db.delete(domain)
                    else:
                        c2 = pickle.loads(data)
                        if name is None:
                            del c2[path]
                        else:
                            del c2[path][name]
                else:
                    raise KeyError("no domain '%s'" % domain)
        finally:
            self._cookies_lock.release()
    def set_cookie(self, cookie):
        db = self._db
        self._cookies_lock.acquire()
        try:
            # store 2-level dict under domain, like {path: {name: value}}
            data = db.get(cookie.domain)
            if data is None:
                c2 = {}
            else:
                c2 = pickle.loads(data)
            if not c2.has_key(cookie.path): c2[cookie.path] = {}
            c3 = c2[cookie.path]
            c3[cookie.name] = cookie
            db.put(cookie.domain, pickle.dumps(c2))
        finally:
            self._cookies_lock.release()
    def __iter__(self):
        return BSDDBIterator(self._db.cursor())
    def _cookies_for_domain(self, domain, request, unverifiable):
        debug("Checking %s for cookies to return", domain)
        if not self.policy.domain_return_ok(domain, request, unverifiable):
            return []

        data = self._db.get(domain)
        if data is None:
            return []
        cookies_by_path = pickle.loads(data)

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
