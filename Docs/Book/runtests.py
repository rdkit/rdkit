#!/usr/bin/env python
"""
  Looks through all ODF and text files in a directory and tries to run
  them as doctests.

  In the ODF files the code samples are recognized as things in
  styles CDT, CDT1, CDTX, C1, or CodeSample

  Modifications to standard doctest behavior :
  1) sample output lines that begin with "#" are ignored. This allows
  output that doctest wouldn't normally "see" (i.e. stuff to stderr or
  things sent to stdout by extension modules) to be shown in the
  documentation without causing tests to fail

  2) everything on sample output lines that follows "#<-" is
  ignored. This allows clarifications in the sample output. NOTE: the
  test script uses the NORMALIZE_WHITESPACE option, so this
  modification will cause everything in the sample output following
  the #<- to be ignored. Here's an example:
    >>> dosomething()
    3.5 7.4   #<- Everything after this point is ignored
    this is all
    ignored as
    well
    >>> 

  Modified by Greg Landrum from the file :
  http://svn.python.org/view/sandbox/trunk/userref/test_examples.py

"""

import os, os.path
import doctest # Custom doctest that checks sys.displayhook
import re
import zipfile
try:
    from xml.etree.cElementTree import iterparse
except ImportError:
    from cElementTree import iterparse    
from cStringIO import StringIO

class ErrorIgnorer(doctest.OutputChecker):
    def check_output(self,want,got,optionflags):
        if want:
            if want[0]=='#' and not got:
                return True
            want = want.split('#<-')[0]
        return doctest.OutputChecker.check_output(self,want,got,optionflags)

def _teststring(data, name="<doctest>", globs=None, verbose=None,
                report=True, optionflags=0, extraglobs=None,
                raise_on_error=False, parser=doctest.DocTestParser()):
    """
    Test examples in the given data string.  Return (#failures, #tests).

    Optional keyword arg "name" gives the name of the test; by default
    uses "<doctest>".

    Optional keyword arg "globs" gives a dict to be used as the globals
    when executing examples; by default, use {}.  A copy of this dict
    is actually used for each docstring, so that each docstring's
    examples start with a clean slate.

    Optional keyword arg "extraglobs" gives a dictionary that should be
    merged into the globals that are used to execute examples.  By
    default, no extra globals are used.

    Optional keyword arg "verbose" prints lots of stuff if true, prints
    only failures if false; by default, it's true iff "-v" is in sys.argv.

    Optional keyword arg "report" prints a summary at the end when true,
    else prints nothing at the end.  In verbose mode, the summary is
    detailed, else very brief (in fact, empty if all tests passed).

    Optional keyword arg "optionflags" or's together module constants,
    and defaults to 0.  Possible values (see the docs for details):

        DONT_ACCEPT_TRUE_FOR_1
        DONT_ACCEPT_BLANKLINE
        NORMALIZE_WHITESPACE
        ELLIPSIS
        IGNORE_EXCEPTION_DETAIL
        REPORT_UDIFF
        REPORT_CDIFF
        REPORT_NDIFF
        REPORT_ONLY_FIRST_FAILURE

    Optional keyword arg "raise_on_error" raises an exception on the
    first unexpected exception or failure. This allows failures to be
    post-mortem debugged.

    Optional keyword arg "parser" specifies a DocTestParser (or
    subclass) that should be used to extract tests from the files.

    Advanced tomfoolery:  teststring runs methods of a local instance of
    class doctest.Tester, then merges the results into (or creates)
    global Tester instance doctest.master.  Methods of doctest.master
    can be called directly too, if you want to do something unusual.
    Passing report=0 to testmod is especially useful then, to delay
    displaying a summary.  Invoke doctest.master.summarize(verbose)
    when you're done fiddling.
    """
    # Assemble the globals.
    if globs is None:
        globs = {}
    else:
        globs = globs.copy()
    if extraglobs is not None:
        globs.update(extraglobs)
    # Make the test runner
    if raise_on_error:
        runner = doctest.DebugRunner(checker=ErrorIgnorer(),
                                     verbose=verbose,
                                     optionflags=optionflags)
    else:
        runner = doctest.DocTestRunner(checker=ErrorIgnorer(),
                                       verbose=verbose,
                                       optionflags=optionflags)
    # Convert the string to a test, and run it.
    test = parser.get_doctest(data, globs, name, None, 0)
    runner.run(test)
    if report:
        runner.summarize()
    return runner.failures, runner.tries

test_options = doctest.ELLIPSIS | doctest.IGNORE_EXCEPTION_DETAIL | doctest.NORMALIZE_WHITESPACE
import __future__
test_globals = dict(absolute_import=__future__.absolute_import,
                    division=__future__.division,
                    with_statement=__future__.with_statement,
                    __name__="__main__")

def _is_text_example(name):
    return re.search("examples\d\d\.txt", name) is not None

def _test_text_file(name):
    return doctest.testfile(name, optionflags=test_options)

def _is_odt_example(name):
    if re.search(".*Python\.odt", name) is not None:
        return True
    elif re.search(".*Book\.odt", name) is not None:
        return True
    return False


_style_attr = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}style-name"
def _in_style(elem, styles):
    if styles is None:
        return True
    try:
        return elem.attrib[_style_attr] in styles
    except KeyError:
        return False

_space_tag = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}s"
_num_spaces_attr = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}c"
_para_tag = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}p"
_span_tag = "{urn:oasis:names:tc:opendocument:xmlns:text:1.0}span"
def _para_text(elem, styles=None):
    if elem.tag != _para_tag:
        return ''
    if not _in_style(elem, styles):
        return '\n'
    para = []
    # print elem
    if elem.text is not None:
        para.append(elem.text)
    for subelem in elem:
        if subelem.tag == _space_tag:
            try:
                num_spaces = int(subelem.attrib[_num_spaces_attr])
            except KeyError:
                pass
            else:
                para.append(" " * num_spaces)
        elif subelem.tag == _span_tag and subelem.text is not None:
            para.append(subelem.text)
        if subelem.tail is not None:
            para.append(subelem.tail)
    if elem.tail is not None:
        para.append(elem.tail)
    if para and not para[-1].endswith('\n'):
        para.append('\n')
    result = ''.join(para)
    if result.rstrip():
         # print result
         return result
    return ''

def _odt_data_to_text(data, styles=None):
    result = []
    for event, elem in iterparse(data):
        result.append(_para_text(elem, styles)) 
    return ''.join(result)

def _get_odt_text(name, styles=None):
    zf = zipfile.ZipFile(name, 'r')
    data = zf.read("content.xml")
    zf.close()
    return _odt_data_to_text(StringIO(data), styles)

_code_styles = "CDT CDT1 CDTX C1 CodeSample".split()
#outF = file('rt.out','w+')
def _test_odt_file(name):
    txt = _get_odt_text(name, _code_styles)
    return _teststring(txt, name=os.path.basename(name),
                            globs=test_globals,
                            optionflags=test_options)

extension_checks = {
    ".txt" : _is_text_example,
    ".odt" : _is_odt_example,
}

extension_tests = {
    ".txt" : _test_text_file,
    ".odt" : _test_odt_file,
}

def abs_walk(walk_dir):
    walk_path = os.path.abspath(walk_dir)
    for path, subdirs, fnames in os.walk(walk_path):
        for fname in fnames:
            abs_name = os.path.join(path, fname)
            __, ext = os.path.splitext(fname)
            yield abs_name, ext

def _test_all(test_dir='.'):
    import sys
    if len(sys.argv)<2:
        for fname, fext in abs_walk(test_dir):
            try:
                is_example = extension_checks[fext]
                test_example = extension_tests[fext]
            except KeyError:
                continue
            if is_example(fname):
                failed, tests = test_example(fname)
                print ("Passed %d of %d tests in %s" %
                       (tests - failed, tests, fname))
    else:
        for fname in sys.argv[1:]:
            __, fext = os.path.splitext(fname)
            try:
                is_example = extension_checks[fext]
                test_example = extension_tests[fext]
            except KeyError:
                continue
            if is_example(fname):
                failed, tests = test_example(fname)
                print ("Passed %d of %d tests in %s" %
                       (tests - failed, tests, fname))

if __name__ == "__main__":
    _test_all()

    
