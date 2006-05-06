#!/usr/bin/env python
#
# $Id: test_happydoc.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $
#
# Time-stamp: <01/02/03 12:29:34 dhellmann>
#
# Copyright Doug Hellmann 2000
#
#                         All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of Doug
# Hellmann not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
# NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

"""Unit test cases for HappyDoc.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: test_happydoc.py,v $',
    'rcs_id'       : '$Id: test_happydoc.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $',
    'creator'      : 'Doug Hellmann <doughellmann@bigfoot.com>',
    'project'      : 'HappyDoc',
    'created'      : 'Sun, 13-Aug-2000 10:16:13 EDT',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 20:59:27 $',
}

#
# Import system modules
#
import sys
import os
import tempfile
import string
from glob import glob

#
# Fix the import path so that we can get to the unit test module.
#
sys.path.append('../../Python/pyunit')

import unittest


#
# Import Local modules
#
import CommandLineApp

#
# Module
#

class HappyDocTest(unittest.TestCase):

    def __init__(self,
                 outputDir='DefaultTestOutputDir',
                 outputSubDir='',
                 methodName='runTest',
                 ):
        unittest.TestCase.__init__(self, methodName)
        self.output_dir = outputDir
        self.output_sub_dir = outputSubDir
        return

    def setUp(self):
        self.happydoc = './happydoc'
        return

    def runHappyDoc(self, testName, modules=(), argStr='-v'):
        #
        # Set up locals
        #
        output_dir = self.output_dir
        output_sub_dir = self.output_sub_dir
        happydoc = self.happydoc
        output_dir = '%(output_dir)s/%(output_sub_dir)s/%(testName)s' % locals()
        modules = string.join(modules, ' ')
        #
        # Build command
        #
        command = '%(happydoc)s -t %(testName)s -d %(output_dir)s %(argStr)s formatter_dateStampFiles=0 %(modules)s 2>&1' % locals()

        print '=' * 80
        print
        print 'COMMAND: %s' % command
        print
        sys.stdout.flush()
        hdRC = os.system(command)
        return hdRC
    
    def checkHelpSyntax(self):
        assert not self.runHappyDoc('checkHelpSyntax', (), argStr='-h'), 'Basic help syntx test failed.'
        return

    def checkHelpManual(self):
        assert not self.runHappyDoc('checkHelpManual', (), argStr='--help'), 'Help manual generation test failed.'
        return

    def checkBasicParser(self):
        assert (not self.runHappyDoc('checkBasicParser',
                                     ('TestCases/test.py',
                                      'TestCases/test_just_docstring.py',
                                      ),
                                     #argStr='-T singlefile',
                                     )
                ), 'Basic parser test failed.'
        return

    def checkDuboisParser(self):
        assert (not self.runHappyDoc('checkDuboisParser',
                                     ('TestCases/mymod.py',),
                                     argStr='-v -p -',
                                     )
                ), 'Paul Dubois submitted parser test failed.'
        return
 
    def checkSingleFileDocset(self):
        assert (not self.runHappyDoc('checkSingleFileDocset',
                                     ('TestCases/test.py',),
                                     argStr='-p - -v -T singlefile')
                ), 'Basic single-file docset test failed.'
        return

    def checkTextStdoutDocset(self):
        assert (not self.runHappyDoc('checkTextStdoutDocset',
                                     ('TestCases/test.py',),
                                     argStr='-v -F text -T stdout')
                ), 'Text to standard-output docset test failed.'
        return

    def checkHTMLStdoutDocset(self):
        assert (not self.runHappyDoc('checkHTMLStdoutDocset',
                                     ('TestCases/test.py',),
                                     argStr='-v -T stdout')
                ), 'HTML to standard-output docset test failed.'
        return

    def checkTextTestModule(self):
        assert (not self.runHappyDoc('checkTextTestModule',
                                     ('TestCases/test.py',),
                                     argStr='-v -F text')
                ), 'Text formatter test failed.'
        return
    
    def checkBasicTwoModules(self):
        assert (not self.runHappyDoc('checkBasicTwoModules',
                                     ('TestCases/test.py', 'CommandLineApp.py')
                                     )
                ), 'Two related modules test failed.'
        return

    def checkBasicImports(self):
        assert (not self.runHappyDoc('checkBasicImports',
                                     ('TestCases/test_import_statements.py',
                                      'CommandLineApp.py',
                                      'prettyast.py',
                                      'parseinfo.py'),
                                     argStr='-v -p ""'
                                     )
                ), 'Imported module id test failed.'
        return

    def checkBasicImportsFromPackages(self):
        assert (not self.runHappyDoc('checkBasicImportsFromPackages',
                                     ('TestCases/test_import_packages',),
                                     argStr='-v -p ""'
                                     )
                ), 'Import from packages test failed.'
        return

    def checkOutputPrefix(self):
        assert (not self.runHappyDoc('checkOutputPrefix',
                                     ('TestCases/test_import_packages',),
                                     argStr='-v -p "" formatter_filenamePrefix="TESTPREFIX_"'
                                     )
                ), 'Formatter output prefix test failed.'
        return

    def checkBasicFunctionParameters(self):
        assert (not self.runHappyDoc('checkBasicFunctionParameters',
                                     ('TestCases/test_function_params.py',
                                      ),
                                     argStr='-v -p ""'
                                     )
                ), 'Function parameter test failed.'
        return
    
    def checkBasicReadme(self):
        assert (not self.runHappyDoc('checkBasicReadme')), 'README test failed.'
        return
    
    def checkSelfDocumentation(self):
        assert (not self.runHappyDoc('checkSelfDocumentation',
                                     ('../HappyDoc',))
                ), 'Full self-documentation test failed.'
        return
    
    def checkSelfDocumentationCompact(self):
        assert (not self.runHappyDoc('checkSelfDocumentationCompact',
                                     ('../HappyDoc',),
                                     argStr='hdformatter_compactHTML="yes"')
                ), 'Full self-documentation with compact output test failed.'
        return
    
    def checkSelfDocumentationDocBook(self):
        assert (not self.runHappyDoc('checkSelfDocumentationDocBook',
                                     ('../HappyDoc',),
                                     argStr='-F sgmldocbook')
                ), 'Full self-documentation in DocBook format output test failed.'
        return
    
    def checkSelfDocumentationDocBookSingleFile(self):
        assert (not self.runHappyDoc('checkSelfDocumentationDocBookSingleFile',
                                     ('../HappyDoc',),
                                     argStr='-F sgmldocbook -T singlefile')
                ), 'Full self-documentation in DocBook format output test failed.'
        return
    
    def checkSelfDocumentationNoDescription(self):
        assert (not self.runHappyDoc('checkSelfDocumentationNoDescription',
                                     ('../HappyDoc',),
                                     argStr='-p -')
                ), 'Full self-documentation test failed.'
        return

    def checkSelfDocumentationAsText(self):
        assert (not self.runHappyDoc('checkSelfDocumentationAsText',
                                     ('../HappyDoc',),
                                     argStr='-v -F text')
                ), 'Full self-documentation-as-text test failed.'
        return

    def checkPrivateNames(self):
        assert (not self.runHappyDoc('checkPrivateNames',
                                     ('TestCases/test_private_names.py',),
                                     argStr='-v --no_private_names -p ""')
                ), 'Documentation without showing private names.'
        return

    def checkIgnoreComments(self):
        assert (not self.runHappyDoc('checkIgnoreComments',
                                     ('TestCases/test_ignore_comments.py',),
                                     argStr='-v --no_comments -p ""')
                ), 'Documentation without showing comments as __doc__ strings.'
        return
    

class OtherWorkingDirTest(HappyDocTest):

    def __init__(self,
                 workingDir='.',
                 outputDir='DefaultTestOutputDir',
                 **nargs
                 ):
        self.dir_stack = None
        self.working_dir = workingDir
        apply(HappyDocTest.__init__, (self,), nargs)
        self.output_dir = os.path.join(os.pardir, 'HappyDoc', outputDir)
        return

    def setUp(self):
        self.happydoc = '../HappyDoc/happydoc'
        return
    
    def runHappyDoc(self, *args, **nargs):
        self.pushDir()
        apply(HappyDocTest.runHappyDoc, (self,) + args, nargs)
        self.popDir()
        return

    def pushDir(self):
        self.dir_stack = (os.getcwd(), self.dir_stack)
        os.chdir(self.working_dir)
        return

    def popDir(self):
        if self.dir_stack:
            top, self.dir_stack = self.dir_stack
            os.chdir(top)
        return

class ZopeTest(OtherWorkingDirTest):
    
    def checkZopeFull(self):
        assert (not self.runHappyDoc('checkZopeFull',
                                     ('../Zope-2-CVS-src',),
                                     argStr='-v')
                ), 'Zope full documentation test failed.'
        return
    
    def checkZopeRoot(self):
        assert (not self.runHappyDoc('checkZopeRoot',
                                     ('../Zope-2-CVS-src',),
                                     argStr='-v -r')
                ), 'Zope full documentation test failed.'
        return

    def checkGadflyParseError(self):
        assert (not self.runHappyDoc('checkGadflyParseError',
                                     ('../Zope-2-CVS-src/lib/python/Products/ZGadflyDA/gadfly/gfdb0.py',),
                                     argStr='-v -r')
                ), 'Gadfly test with parse-error failed.'
        return

    def checkZEOParseError(self):
        assert (not self.runHappyDoc('checkZEOParseError',
                                     ('../Zope-2-CVS-src/lib/python/ZEO/zrpc.py',),
                                     argStr='-v -r')
                ), 'ZEO test with parse-error failed.'
        return
        
        
        
def ZopeTestFactory(**nargs):
    nargs['workingDir'] = nargs.get('workingDir', '../Zope-2-CVS-src')
    return apply(ZopeTest, (), nargs)

class HappyDocBugRegressionTest(HappyDocTest):

    def __init__(self,
                 outputDir='DefaultTestOutputDir',
                 outputSubDir='',
                 methodName='',
                 ):
        HappyDocTest.__init__(self,
                              outputDir=outputDir,
                              outputSubDir=outputSubDir,
                              methodName='checkBugReport%s' % methodName)
        return

    def checkBugReport(self, bugId):
        assert not self.runHappyDoc(bugId, ('TestCases/test_bug%s.py' % bugId,
                                            ),
                                    argStr='-p -'), \
                                    'Check for bug %s failed.' % bugId

    def __getattr__(self, name):
        if name[:14] == 'checkBugReport':
            id = name[14:]
            return lambda bug=id, s=self: s.checkBugReport(bug)
        raise AttributeError(name)
    
    
class TestCaseDriver(CommandLineApp.CommandLineApp):
    "Drive the test cases for HappyDoc."

    def optionHandler_t(self, testSetName):
        "Specify a testSetName to run."
        self.test_set = testSetName
        return

    def optionHandler_d(self, outputDir):
        "Specify an output directory for the test output files."
        self.output_dir = outputDir
        return

    def appInit(self):
        self.optionHandler_t('default')
        self.optionHandler_d('../HappyDocRegressionTest/SimpleTestOutput')
        return

    def buildTests(self):
        "Create the test suites."

        #
        # Test suite to run all tests
        #
        self.allTestSuite = unittest.TestSuite()
        add_all_tests = self.allTestSuite.addTest

        bug_ids = map(lambda x:x[18:-3], glob('TestCases/test_bug*.py'))
        bug_ids.sort()
        
        for suite_name, suite_type, definitions in (

            #
            # The default test suite
            #
            ( 'default', HappyDocTest, (
            'checkBasicImportsFromPackages',
            ) ),
            
            #
            # Test help generation
            #
            ( 'help', HappyDocTest, ( 'checkHelpSyntax',
                                      'checkHelpManual',
                                      ) ),

            #
            # Test basic pieces of the parser and formatters
            #
            ( 'basic', HappyDocTest, ( 'checkBasicParser',
                                       'checkBasicReadme',
                                       'checkBasicTwoModules',
                                       'checkBasicImports',
                                       'checkBasicFunctionParameters',
                                       'checkPrivateNames',
                                       'checkIgnoreComments',
                                       ) ),
              
            #
            # Test processing the HappyDoc code and generating docs
            #
            ( 'self', HappyDocTest, ( 'checkSelfDocumentation',
                                      'checkSelfDocumentationNoDescription',
                                      'checkSelfDocumentationAsText',
                                      'checkSelfDocumentationCompact',
                                      ) ),

            #
            # Test the formatter_textfile.py module.
            #
            ( 'text', HappyDocTest, ( 'checkTextTestModule',
                                      'checkSelfDocumentationAsText',
                                      ) ),
            
            #
            # Test the various docset formats
            #
            ( 'docset', HappyDocTest, ( 'checkSingleFileDocset',
                                        'checkTextStdoutDocset',
                                        'checkHTMLStdoutDocset',
                                        ) ),
            
            #
            # Test against the Zope source code
            #
            ( 'zope_root', ZopeTestFactory, ('checkZopeRoot',) ),
            
            ( 'zope', ZopeTestFactory, ('checkZopeFull',) ),

            ( 'ZEO', ZopeTestFactory, ('checkZEOParseError',
                                        ) ),

            #
            # Test against submitted parser error code.
            #
            ( 'parse_error', HappyDocTest, ('checkDuboisParser',),
              ),
            
            #
            # Test against the Gadfly code which has
            # parse errors.
            #
            ( 'zope_parse_error', ZopeTestFactory, ( 'checkGadflyParseError',
                                                ),
              ),

            #
            # Check tests related to bug reports
            #
            ( 'bugs', HappyDocBugRegressionTest, bug_ids,
              ),

            #
            # The hdformatter filename prefix test suite
            #
            ( 'output_prefix', HappyDocTest, (
            'checkOutputPrefix',
            ) ),

            #
            # Packages
            #
            ( 'package', HappyDocTest, (
            'checkBasicImportsFromPackages',
            ) ),

            #
            # DocBook
            #
            ( 'docbook', HappyDocTest, (
            'checkSelfDocumentationDocBook',
            'checkSelfDocumentationDocBookSingleFile',
            ) ),
            

            ):
            test_suite = unittest.TestSuite()
            setattr( self, '%sTestSuite' % suite_name, test_suite )
            for test_name in definitions:
                new_test = suite_type(outputDir=self.output_dir,
                                      outputSubDir=suite_name,
                                      methodName=test_name,
                                      )
                test_suite.addTest( new_test )
            add_all_tests( test_suite )

        return

    def main(self, *args):
        "Run the required test suites."

        if args:
            raise ValueError('Unhandled arguments!', args)

        self.buildTests()

        actual_test_suite = unittest.TestSuite()

        try:
            test_suite = getattr(self, '%sTestSuite' % self.test_set)
        except AttributeError:
            self.errorMessage('Unrecognized test suite "%s" enabled.' % self.test_set)
            self.errorMessage('Skipping.')
        else:
            actual_test_suite.addTest(test_suite)
                    
        #
        # Run the test suite
        #

        print '=' * 80
        print 'START'
        print '-' * 80
        print

        runner = unittest.TextTestRunner(sys.stdout)
        runner.run(actual_test_suite)

        print
        print '-' * 80
        print 'FINISH'
        print '=' * 80

        return
    

if __name__ == '__main__':
    try:
        TestCaseDriver().run()
    except TestCaseDriver.HelpRequested:
        pass
    
