#!/usr/bin/env python
#
# Time-stamp: <01/02/03 12:55:20 dhellmann>
#
# COPYRIGHT
#
#   Permission to use, copy, modify, and distribute this software and
#   its documentation for any purpose and without fee is hereby
#   granted, provided that the above copyright notice appear in all
#   copies and that both that copyright notice and this permission
#   notice appear in supporting documentation, and that the name of Doug
#   Hellmann not be used in advertising or publicity pertaining to
#   distribution of the software without specific, written prior
#   permission.
# 
# DISCLAIMER
#
#   DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#   INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
#   NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY SPECIAL, INDIRECT OR
#   CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
#   OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
#   NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
#   CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 


"""A base class for file-based formatters for HappyDoc.

$Id: hdformatter_file.py,v 1.1 2001/03/21 21:58:54 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: hdformatter_file.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 17:56:22 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 21:58:54 $',
    }

#
# Import system modules
#
import os
import string
import types

#
# Import Local modules
#
import happyformatter
from hdpathmanagement import rmkdir, applyPrefixToPath

#
# Module
#

    
class FileBasedFormatter(happyformatter.HappyFormatterBase):
    "Output documentation to file(s)."

    def __init__(self, docSet, filenamePrefix='', **extraNamedParameters):
        """Initialize a FileBasedFormatter.

        Parameters

          docSet -- The documentation set controlling the formatter.
          
          filenamePrefix -- A prefix to append to the base names of
          files and directories being created.  This is useful for
          situations where the names which would be automatically
          generated might cause a name clash or conflict.

          extraNamedParameters -- Parameters specified by name which
          were not supported by a subclass initialization.
          
        """
        #
        # Store parameters
        #
        self._filename_prefix = filenamePrefix
        #
        # Initialize the base class
        #
        apply(happyformatter.HappyFormatterBase.__init__,
              (self, docSet,),
              extraNamedParameters)
        return

    def openOutput(self, name, title1, title2=None):
        "Open the named output destination and give the title."
        rmkdir(os.path.dirname(name))
        f = open(name, 'wt')
        if not hasattr(self, 'open_root_file'):
            self.open_root_file = f
        return f

    def closeOutput(self, output):
        "Close the output handle."
        output.close()
        return
        
    def getOutputNameForFile(self, filename, usePrefix=0):
        """
        Return the base name of the thing to which output should be
        written for a file.  This is usually a file name, but could be
        anything understood by the formatter as a name.  If infoObject
        is None, return the name for a root node for this formatter.
        """
        #
        # Remove preceding slashes to make name relative
        #
        while filename[0] in './':
            filename = filename[1:]
        #
        # Apply the path prefix, if required
        #
        if usePrefix:
            filename = applyPrefixToPath(filename, self._filename_prefix)
        #
        # Set the correct extension for the output file
        #
        extension = self.getFilenameExtension()
        #
        # Build the name
        #
        basename, ext = os.path.splitext(filename)
        if ext == '.py':
            name = '%s.%s' % (filename, extension)
        else:
            name = '%s.%s' % (basename, extension)
        return name
    
    def getOutputNameForObject(self, infoObject):
        """
        Return the base name of the thing to which output should be written
        for an info source.  This is usually a file name, but could
        be anything understood by the formatter as a name.  If
        infoObject is None, return the name for a root node for this
        formatter.
        """
        #print 'ZZZ getOutputNameForObject(%s)' % infoObject
        if type(infoObject) == types.StringType:
            #print 'ZZZ string'
            name = infoObject
        elif type(infoObject) == types.FileType:
            #print 'ZZZ file'
            name = infoObject.name
        elif infoObject is not None:
            #print 'ZZZ infoObject'
            name = self.getOutputNameForFile(
                infoObject.getFullyQualifiedName(self._filename_prefix)
                )
            #print '\n\nFILE for %s \nis %s\n\n' % (infoObject.getName(),
            #                                       name)
        else:
            name = self.getRootNodeName()
        #print 'ZZZ =>%s' % name
        return name


    def getLocalOutputNameForObject(self, infoObject):
        """
        Return a local reference to base name of the thing to which
        output should be written for an info source.  This is usually
        a file name, but could be anything understood by the formatter
        as a name.  If infoObject is None, return the name for a root
        node for this formatter.
        """
        extension = self.getFilenameExtension()
        if infoObject is not None:
            name = '%s.%s' % ( infoObject.getQualifiedName(self._filename_prefix),
                               extension )
        else:
            name = self.getRootNodeName()
        return name


    def getFullOutputNameForObject(self, infoObject):
        "Get the full name, including path, to the object being output."
        #print 'AAA getFullOutputNameForObject(%s)' % infoObject
        if self._docset._base_directory:
            name = os.path.join( self._docset._base_directory,
                                 self.getOutputNameForObject(infoObject)
                                 )
        else:
            name = self.getOutputNameForObject(infoObject)
        #print 'AAA =>%s' % name
        return name

    def getFullOutputNameForFile(self, filename, usePrefix=0):
        "Get the full name, including path, to the filename to convert."
        if self._docset._base_directory:
            return os.path.join( self._docset._base_directory,
                                 self.getOutputNameForFile(filename, usePrefix)
                                 )
        else:
            return self.getOutputNameForFile(filename)


    def getRootNodeName(self):
        "Return the name of the root node for the documentation tree."
        self._requiredOfSubclass('getRootNodeName')
        return
