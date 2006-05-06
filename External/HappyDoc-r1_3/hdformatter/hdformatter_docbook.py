#!/usr/bin/env python
#
# $Id: hdformatter_docbook.py,v 1.1 2001/03/21 21:58:54 glandrum Exp $
#
# Time-stamp: <01/02/03 13:07:46 dhellmann>
#
# Copyright 2001 Balazs Scheidler
#
#
#                         All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of Balazs
# Scheidler not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# BALAZS SCHEIDLER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
# NO EVENT SHALL BALAZS SCHEIDLER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

"""Output documentation in SGML using DocBook markup

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: hdformatter_docbook.py,v $',
    'rcs_id'       : '$Id: hdformatter_docbook.py,v 1.1 2001/03/21 21:58:54 glandrum Exp $',
    'creator'      : 'Balazs Scheidler <bazsi@balabit.hu>',
    'project'      : 'HappyDoc',
    'created'      : 'Sat, 03-Feb-2001 12:53:37 EST',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 21:58:54 $',
}

#
# Import Local modules
#
import hdformatter_file
import StructuredText
import indentstring


#
# Module
#

def entryPoint():
    "Return information about this module to the dynamic loader."
    return {
        'name':'sgmldocbook',
        'factory':DocBookFormatter,
        }


class DocBookFormatter(hdformatter_file.FileBasedFormatter):
    def __init__(self, docSet, **configuration):
        self._section_level_counter = 1
        self.debug = 1
        apply(hdformatter_file.FileBasedFormatter.__init__, (self, docSet),
              configuration)
        return
    
    def openOutput(self, name, title1, title2=""):
        """Write the formatting for a file header to the open file."""
        f = hdformatter_file.FileBasedFormatter.openOutput(self,
                                                           name,
                                                           title1)
        self.fileHeader(title1, title2, f)
        return f

    def closeOutput(self, output):
        "Close the 'output' handle."
        self.fileFooter(output)
        output.close()
        return

    def fileHeader(self, title1, title2, output):
        self.comment('file_header', output)
        return

    def fileFooter(self, output):
        self.comment('file_footer', output)
        return

    # string handling

    def writeText(self, text, output, quote=0):
        text = self._unquoteString(text)
        #html = str(StructuredText.HTML(text))
        self.writeRaw(text, output)

    def formatCode(self, text):
        "Format 'text' as source code and return the new string."
        formatted_text = '<programlisting>\n%s\n</programlisting>\n' % text
        return formatted_text

    def formatKeyword(self, text):
        "Format 'text' as a keyword and return the new string."
        formatted_text = '<literal>%s</literal>' % text
        return formatted_text

    # structure handling


    # simple lists
    def listHeader(self, output, title=None, allowMultiColumn=1):
        """Output 'title' as a heading for a list.  If 'allowMultiColumn' is
        true, set up the list to have more than one column.
        """
        self.writeRaw('<formalpara>\n<title>%s</title>\n<para><itemizedlist>\n' % title, output)
        return

    def listItem(self, output, text):
        "Format and output the 'text' as a list element."
        self.writeRaw('<listitem><para>%s</para></listitem>\n' % text, output)
        return

    def listFooter(self, output):
        "Write the closing footer for a list to the 'output'."
        self.writeRaw('\n</itemizedlist></para></formalpara>', output)
        return

    # descriptive lists
    def descriptiveListHeader(self, output, title):
        "Write the 'title' as the heading for a descriptive list to the 'output'."
        self.writeRaw('<formalpara><title>%s</title>\n<para><variablelist>\n' \
                      % title, output)
        self.comment('descriptive list header', output)
        return

    def descriptiveListItem(self, output, item, description):
        "Format and write the 'item' and 'description' for a descriptive list to the 'output'."
        self.writeRaw('<varlistitem><term>%s</term><listitem><para>' % item,
                      output)
        self.writeText(description, output)
        self.writeRaw('</para></listitem></varlistitem>\n', output)
        return

    def descriptiveListFooter(self, output):
        "Write the closing footer for a descriptive list to the 'output'."
        self.writeRaw('</variablelist></para></formalpara>\n', output)
        return

    # headers

    def sectionHeader(self, output, title):
        "Write a general purpose section openning title to the 'output'."
        self.writeRaw('<sect%d>\n<title>%s</title>' \
                      % (self._section_level_counter, title),
                      output)
        return
       
    def sectionFooter(self, output):
        "Write a general purpose section closing footer to the 'output'."
        self.writeRaw('</sect%d>' % self._section_level_counter, output)
        return

    def itemHeader(self, output, infoObject):
        "Write a section openning header for an 'infoObject' to the 'output'."
        name = infoObject.getName()
        self.sectionHeader(output, name)
        return

    def itemFooter(self, output):
        "Write a section closing footer to the 'output'."
        self.sectionFooter(output)
        return

    def pushSectionLevel(self, output):
        self._section_level_counter = self._section_level_counter + 1
        return

    def popSectionLevel(self, output):
        self._section_level_counter = self._section_level_counter - 1
        return

    # misc

    def dividingLine(self, output, fill='-'):
        "Write a sectional dividing line made up of repeated 'fill' characters to the 'output'."
        #output.write('<hr>\n')
        return

    def comment(self, text, output):
        """Output text as a comment."""
        if self.debug: self.writeRaw('<!-- %s -->\n' % text, output)
        return

    def indent(self, output):
        self.comment('indent', output)
        return

    def dedent(self, output):
        self.comment('dedent', output)
        return


    # refererences

    def getReference(self, infoSource, relativeSource):
        """Returns a reference to the 'infoSource' from 'relativeSource'.
        """
        #print '\ngetReference(', infoSource, ',', relativeSource, ')'
        #link = computeRelativeLink(relativeSource,
        #                           self.getOutputNameForObject(infoSource)
        #                           )
        info = {
            'name':self.getNameForInfoSource( infoSource ),
            }
        ref = '<xref linkend="%(name)s">' % info
        return ref

    def getNamedReference(self, infoSource, name, relativeSource):
        info = {
            'name':infoSource.getName(),
            }
        ref = '<xref linkend="%(name)s">' % info
        return ref

    def getInternalReference(self, infoSource):
        """Returns a reference to 'infoSource' within the current document.
        """
        info = {
            'name':infoSource.getName(),
            }
        ref = '<xref linkend="%(name)s">' % info
        return ref
   
    def getPythonReference(self, moduleName):
        """Returns a reference to 'moduleName' documentation on the
       
        "Python.org":http://www.python.org documentation site.
        """
        if moduleName in self.sys_modules:
            return '<ulink url="http://www.python.org/doc/current/lib/module-%(moduleName)s.html">%(moduleName)s</ulink>' % locals()
        else:
            return moduleName

    def getFilenameExtension(self):
        "Returns the extension for creating output files."
        return 'sgml'

    def getRootNodeName(self):
        return 'book.sgml'
