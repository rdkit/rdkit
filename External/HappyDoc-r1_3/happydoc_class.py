#!/usr/bin/env python
#
# $Id: happydoc_class.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $
#
# Time-stamp: <01/02/03 11:17:44 dhellmann>
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

"""HappyDoc Class for HappyDoc Application.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: happydoc_class.py,v $',
    'rcs_id'       : '$Id: happydoc_class.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $',
    'creator'      : 'Doug Hellmann <doughellmann@bigfoot.com>',
    'project'      : 'HappyDoc',
    'created'      : 'Sun, 13-Aug-2000 11:27:00 EDT',

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
import pprint
import os
import glob
import sys
import types
import string
import re
import parser
import traceback

#
# Import Local modules
#
from parseinfo import getDocs
import parseinfo
from CommandLineApp import CommandLineApp
from prettyast import astListFixNames
import docset
import hdformatter
from cvsversion import cvs_product_version
from hdpathmanagement import rmkdir
import optiontools

#
# Module
#

True = 1
False = None
    
class HappyDoc(CommandLineApp):
    """
    HappyDoc is a documentation generation/extraction tool which does
    not depend on being able to import modules.

    This app is based on the Demos/parser/example.py module distributed
    with the Python source distribution.

    """

    include_private_names = True
    include_comments = True
    output_directory = './doc'
    output = None
    author_name = 'Doug Hellmann <doughellmann@bigfoot.com>'
    app_home = 'http://HappyDoc.sourceforge.net/'
    package_description_file = 'README.txt'
    recurse_into_subdirs=True

    docset_type = None
    docset_title = 'HappyDoc Generated Documentation'

    #
    # Define the formats supported
    #
    supported_formats = hdformatter.plugins

    #
    # Define the documentation set types supported
    #
    supported_docset_types = docset.plugins
    
    def appInit(self):
        self._app_name = self.__class__.__name__
        self._app_version = cvs_product_version()
        
        self.statusMessage('%s version %s' % (self._app_name,
                                              self._app_version))
        
        self.set_docset_type('docset')
        self.set_format('htmltable')
        
        self._ignore_dirs = []
        self.add_ignore_directory('CVS', 'dist', 'build')
        return

    def add_ignore_directory(self, *dirNames):
        "Add one or more directory names to the list which should be ignored."
        for dir_name in dirNames:
            if dir_name not in self._ignore_dirs:
                self._ignore_dirs.append(dir_name)
        pattern_string = r'^.*/(%s)$' % string.join(self._ignore_dirs,
                                                    '|')
        self._ignore_dir_name = re.compile(pattern_string).match
        return
    
    def set_format(self, format):
        "Set the formatter to be used."
        self.format = format
        try:
            self.formatter = self.supported_formats[format]
        except KeyError:
            raise ValueError('format must be one of %s' \
                             % self.supported_formats.keys(),
                             format)
        return

    def set_docset_type(self, docset_type):
        "Set the docset to be used."
        self.docset_type = docset_type
        try:
            self.docset = self.supported_docset_types[docset_type]
        except KeyError:
            raise ValueError('docset_type must be one of %s' % \
                             self.supported_docset_types.keys(),
                             docset_type)
        return

    def optionHandler_no_private_names(self):
        "Do not include names beginning with _."
        self.include_private_names = False
        return

    def optionHandler_no_comments(self):
        """Do not include comment text as though it was
           a __doc__ string.
        """
        self.include_comments = False
        return
    
    def optionHandler_d(self, outputDirectory):
        "Specify an outputDirectory."
        self.output_directory = outputDirectory
        return

    def optionHandler_r(self):
        "Disable recursion into subdirectories."
        self.recurse_into_subdirs = False
        return

    def optionHandler_o(self):
        "Specify that output should go to stdout."
        self.set_docset_type('stdout')
        return

    def optionHandler_t(self, title):
        "Specify a title for the documentation set."
        self.docset_title = title
        return

    def optionHandler_p(self, packageDescriptionFile):
        """Specify a file with a description of the package.

        The default packageDescriptionFile is README.txt.
        """
        self.package_description_file = packageDescriptionFile
        return

    def optionHandler_author(self, authorNameAndEmail):
        """Specify the author identification to be inserted for
        references.
        """
        self.author_name = authorNameAndEmail
        return

    def optionHandler_F(self, format):
        "Specify the output format."
        self.set_format(format)
        return

    def optionHandler_T(self, docset_type):
        "Specify the documentation set type."
        self.set_docset_type(docset_type)
        return

    def showVerboseSyntaxHelp(self):
        "Overloaded to show supported docset and format types."
        CommandLineApp.showVerboseSyntaxHelp(self)

        print 'SUPPORTED FORMATS for -F Option:\n'
        items = self.supported_formats.items()
        items.sort()
        for i in items:
            print '  FORMATTER TYPE %s: %s' % (i[0], i[1].__doc__)

        print '\nSUPPORTED DOCSET TYPES for -T Option:\n'
        items = self.supported_docset_types.items()
        items.sort()
        for i in items:
            print '  DOCSET TYPE: %s: %s' % (i[0], i[1].__doc__)

        print
        return

    
    def main(self, *args):
        #
        # Debug info about where the docsets and hdformatters come from
        #
        self.statusMessage('Docsets list from %s' % docset.__path__[0], 1)
        self.statusMessage('Formatters from %s' % hdformatter.__path__[0], 1)
        #
        # Find DocSet arguments
        #
        self.statusMessage('Looking for docset parameters', 2)
        args, docset_params = optiontools.getParameters('docset', args)
        self.statusMessage('DEBUG: Docset parameters:', 3)
        for p, v in docset_params.items():
            self.statusMessage('DEBUG: \t%s:%s' % (p,v), 3)
        #
        # Find Formatter parameters
        #
        self.statusMessage('Looking for formatter parameters', 2)
        args, formatter_params = optiontools.getParameters('formatter', args)
        formatter_params['appVersion'] = formatter_params.get('appVersion',
                                                              self._app_version)
        self.statusMessage('DEBUG: Formatter parameters:', 3)
        for p, v in formatter_params.items():
            self.statusMessage('DEBUG: \t%s:%s' % (p,v), 3)
        #
        # Get the list of modules to input
        #
        input_modules = args
        #
        # Create output directory
        #
        if not self.output:
            od = self.output_directory
            self.statusMessage('Output directory is %s' % self.output_directory, 2)
            if (od[0] == '.'):
                od = os.path.join( os.getcwd(), od )
                self.statusMessage('Setting output directory to %s' % od, 2)
            od = os.path.normpath(od)
            self.statusMessage('Creating output directory %s' % od, 2)
            rmkdir(od)
            self.output_directory = od
        #
        # Create the docset
        #
        docset_init_params = {
            'formatterFactory':self.formatter,
            'parserFunc':parseinfo.getDocs,
            'inputModuleNames':input_modules,
            
            'appHome':self.app_home,
            'appName':self._app_name,
            'appVersion':self._app_version,
            'author':self.author_name,
            'baseDirectory':self.output_directory,
            'descriptionFilename':self.package_description_file,
            'formatterParameters':formatter_params,
            'ignoreDirFunc':self._ignore_dir_name,
            'includeComments':self.include_comments,
            'includePrivateNames':self.include_private_names,
            'statusMessageFunc':self.statusMessage,
            'title':self.docset_title,
            'useRecursion':self.recurse_into_subdirs,

            }
        docset_init_params.update(docset_params)
        parsed_modules = apply( self.docset, (), docset_init_params)
        #
        # Tell the docset to output its results
        #
        parsed_modules.write()
        return

