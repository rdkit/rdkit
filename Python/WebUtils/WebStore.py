#
# copyright (c) 2003 Rational Discovery LLC
#  All Rights Reserved
#
""" A simple, primitive, shelve-based data storage system for persistence
  on web pages.

"""
import sys,cgi,os,time,exceptions
import Cookie
from utils import GUIDs
class TrappedError(exceptions.Exception):
  """ exception class for trapping web errors

  """
  def __init__(self,args=None):
    self.args=args

cookieName = 'RD-Online-Client-ID'
_version = "0.1.0"
cmpdsPerPage = 10

def GetIDCookie():
  """ returns the value of the current ID cookie, _None_ if it's not set

  """
  c = Cookie.SmartCookie()
  try:
    envCookie = os.environ['HTTP_COOKIE']
  except:
    return None

  c.load(envCookie)
  if c.has_key(cookieName):
    v = c[cookieName].value
    return v
  else:
    return None

def ClearIDCookie(req):
  """ clears the current ID Cookie and updates the headers

  """
  c = Cookie.SmartCookie()
  c[cookieName]=None
  if req is not None:
    cOut = c.output(header='',sep='')
    try:
      h = req.headers_out['Set-Cookie']
    except KeyError:
      req.headers_out['Set-Cookie'] = cOut 
    else:
      req.headers_out['Set-Cookie'] = '%s %s'%(h,cOut)
  else:
    print c

def SetIDCookie(req):
  """ sets the value of the ID cookie

  """
  c = Cookie.SmartCookie()
  c[cookieName]=GUIDs.getGUID()
  if req is not None:
    cOut = c.output(header='',sep='')
    try:
      h = req.headers_out['Set-Cookie']
    except KeyError:
      req.headers_out['Set-Cookie'] = cOut 
    else:
      req.headers_out['Set-Cookie'] = '%s %s'%(h,cOut)
  else:
    print c

  return c[cookieName].value

defaultHtmlHeader="""<html>
<head>
  <link rel="stylesheet" type="text/css" href="/RD.css">
  <title>%(title)s</title>
  %(extras)s
</head>
<body>
  <table border=0>
  <tr>
   <td><IMG SRC="/images/RD-Logo.very-sm.jpg" ALT="RD Logo"></td>
   <td><h1>%(title)s</h1></td>
  </tr>
  </table>
   <HR NOSHADE>
"""
def ConstructHtmlHeader(title='',extraHeadBits=''):
  """ returns an HTML header for our templated pages

  """
  argD = {'title':title,'extras':extraHeadBits}
  res = defaultHtmlHeader%argD
  return res

defaultHtmlFooter="""
%(restart)s
%(logout)s
<p><hr><table border=0 width="100%%"><tr>
  <td align=left><i>Version:</i> %(version)s</td>
  <td align =center>Copyright (C) 2002 Rational Discovery LLC</td>
  <td align=right><IMG SRC="/images/RDOnline-Logo2.jpg" ALT="RDOnline Logo"></td>
  </tr></table>
  </body></html>
"""
defaultLogoutText="""<a href="Logout.py">Logout</a>"""
defaultRestartText="""<a href="Init.py">Start Over</a>"""
def ConstructHtmlFooter(includeRestart=1,logoutText=None,restartText=None):
  """ returns an HTML footer for our templated pages

  """
  if logoutText is None:
    logoutText=defaultLogoutText
  if restartText is None:
    restartText=defaultRestartText
    
  argD={'version':str(_version),'restart':'','logout':logoutText}
  if includeRestart:
    argD['restart'] = restartText
  res = defaultHtmlFooter%argD
  return res

def ConstructCookieErrorPage(startPage="Init.py"):
  """ returns the text for a cookie error page

    *(duh)*

  """
  page = ConstructHtmlHeader('Error')
  args = { 'startPage':startPage }
  msg = """
<h3>Cookie Error</h3>
<p>The tool did not find an appropriate cookie set in your browser.
This could be due to a number of reasons:<ol>
<li>You did not start at the <a href="%(startPage)s">start page</a>
<li>You have cookies disabled in your browser.
<li>You have logged out and tried to resume a previous session.
</ol>
The tool needs to be able to set a cookie in order to be able to keep track
of you whilst you use it.  This cookie is set when you visit the
<a href="%(startPage)s">start page</a>, which is where you <b>must</b>
start each session.
<p><hr><center><b><a href="%(startPage)s">Go to the start page</a></b></center>
"""%args
  page = page + msg
  page = page + ConstructHtmlFooter(includeRestart=0)
  return page

def PostHtml(req,page):
  """ posts (by printing) the html passed in.  also sets the content type

  """
  if req is not None:
    req.content_type = 'text/html'
    req.send_http_header()
    req.write(page)
  else:
    print 'Content-type: text/html'
    print ''
    print page

def PostError(msg):
  """ creates an error page

  """
  page = ConstructHtmlHeader('Error')
  page = page + str(msg)
  page = page + ConstructHtmlFooter(includeRestart=0)
  PostHtml(None,page)

def FormFromReq(req):
  """ returns the _cgi.FieldStorage()_

  """
  return cgi.FieldStorage()

def StatusOutput(line):
  """ used to do status output when updating web pages...

    output is dumped to _sys.stdout_ and then flushed to attempt to disable
    buffering.

  """
  sys.stdout.write(line)
  sys.stdout.flush()
