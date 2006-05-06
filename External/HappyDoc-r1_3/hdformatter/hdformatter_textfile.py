#!/usr/bin/env python
#
# Time-stamp: <01/02/03 12:55:39 dhellmann>
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


"""Output documentation in plain ASCII text format.

$Id: hdformatter_textfile.py,v 1.1 2001/03/21 21:58:54 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: hdformatter_textfile.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 17:58:02 EDT',
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
import string
import os
import types

#
# Import Local modules
#
import hdformatter_file

#
# Module
#

def entryPoint():
    "Return information about this module for the dynamic loader."
    return {
        'name':'text',
        'factory':TextFormatter,
        }



class TextFormatter(hdformatter_file.FileBasedFormatter):
    """Output documentation in plain ASCII text format.

    Parameters
          
       filenamePrefix -- A prefix to append to the base names of files
                         and directories being created.  This is
                         useful for situations where the names which
                         would be automatically generated might cause
                         a name clash or conflict.
      
    """

    def writeText(self, text, output, quote=1):
        'Output the text string to the information source.'
        if not text:
            return
        text = string.strip(text)
        if text:
            lines = string.split(text, '\n')
            for line in lines:
                line = string.strip(line)
                if line and line[0] == ' ':
                    line = line[8:]
                elif line and line[0] == '\t':
                    line = line[1:]
                self.writeRaw('\t%s\n' % line, output)
            self.writeRaw('\n', output)
        return

    def writeExceptionListForFunction(self, output, function, listHeader):
        """Write the list of exceptions raised by a function.

        Parameters

          output -- Where to write.

          function -- FunctionInfo from parseinfo module.

          listHeader -- Header for list being generated.

        """
        self.listHeader(output, listHeader)
        exception_names = function.getExceptionNames()
        exception_names.sort()
        output_reduced_name = output.name[len(self._docset.getBaseDirectory())+1:]
        for name in exception_names:
            exception_class = self._docset.getClassInfo(name)
            if exception_class:
                ref = self.getReference(exception_class,
                                        output_reduced_name)
            else:
                ref = self.getPythonReference( name )
            self.listItem(output, self.formatCode(ref))
        self.listFooter(output)
        return
    
    ##
    ## OUTPUT DESTINATIONS
    ##

    def getFilenameExtension(self):
        'Returns the extension to use when creating files for this formatter.'
        return 'txt'

    def getRootNodeName(self):
        "Return the name of the root node for the documentation tree."
        return 'toc.txt'
    
    ##
    ## STRUCTURED OUTPUT METHODS
    ##

    def listHeader(self, output, title=None, allowMultiColumn=0):
        if title:
            self.writeRaw('\n-- %s --\n\n' % title, output)
        #else:
        #    self.writeRaw('\n-- ** --\n\n', output)
        return

    def listItem(self, output, text):
        self.writeRaw('\t* %s\n' % text, output)
        return

    def listHeader(self, output, title=None, allowMultiColumn=1):
        if title:
            self.writeRaw('==> %s <==\n' % title, output)
        return

    def listFooter(self, output):
        self.writeRaw('\n', output)
        return

    def descriptiveListHeader(self, output, title):
        self.listHeader(output, title)
        return

    def descriptiveListItem(self, output, item, description):
        self.writeRaw('\t%s : %s\n\n' % (item, string.strip(description)), output)
        return

    def descriptiveListFooter(self, output):
        self.listFooter(output)
        return

    def sectionHeader(self, output, title):
        self.writeRaw('\n** %s **\n\n' % title, output)
        return

    def sectionFooter(self, output):
        self.writeRaw('\n', output)
        return

    def itemHeader(self, output, infoObject):
        title = infoObject.getName()
        self.writeRaw('\n== %s ==\n\n' % title, output)
        return

    def itemFooter(self, output):
        self.writeRaw('\n', output)
        return

    def dividingLine(self, output, fill='-', span=80):
        '''Output a sectional dividing line.

        Parameters:

            output -- destination for written output

            fill="-" -- character to use to draw the line

            span=80 -- width of line to draw
            
        '''
        span = span / len(fill)
        self.writeRaw(fill * span, output)
        self.writeRaw('\n\n', output)
        return
    
    def pushSectionLevel(self, output):
        pass

    def popSectionLevel(self, output):
        pass

    def getPythonReference(self, moduleName):
        if moduleName in self.sys_modules:
            return moduleName + ' (standard module)'
            #return 'http://www.python.org/doc/current/lib/module-%(moduleName)s.html' % locals()
        else:
            return moduleName
    
    def getReference(self, infoSource, relativeSource):
        #
        # Figure out the name of the infoSource
        #
        if type(infoSource) == types.FileType:
            name = os.path.basename(infoSource.name)
        elif type(infoSource) == types.StringType:
            name = os.path.basename(infoSource)
        else:
            name = infoSource.getName()
            
        info = {
            'name':name,
            'link':self.getOutputNameForObject(infoSource),
            }
        ref = '%(name)s (see %(link)s)' % info
        return ref

    def getLocalReference(self, infoSource):
        info = {
            'name':infoSource.getName(),
            }
        ref = '%(name)s (see %(name)s)' % info
        return ref

    def getInternalReference(self, infoSource):
        info = {
            'name':infoSource.getName(),
            }
        ref = '%(name)s (see %(name)s)' % info
        return ref

    def getNamedReference(self, infoSource, name, relativeSource):
        info = {
            'name':infoSource.getName(),
            'target':name,
            }
        ref = '%(target)s (see %(target)s in %(name)s)' % info
        return ref
        

    def indent(self, output):
        return

    def dedent(self, output):
        return
    
    
    def formatCode(self, text):
        return text

    def formatKeyword(self, text):
        return text
    
    def comment(self, text, output):
        return
    
