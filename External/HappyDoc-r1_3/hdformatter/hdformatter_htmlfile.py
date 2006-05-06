#!/usr/bin/env python
#
# Time-stamp: <01/02/03 13:52:32 dhellmann>
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


"""Output documentation using HTML with tables.

$Id: hdformatter_htmlfile.py,v 1.2 2001/03/21 22:29:33 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: hdformatter_htmlfile.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 17:58:48 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.2 $',
    'date':'$Date: 2001/03/21 22:29:33 $',
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
import hdformatter_file
import StructuredText
import indentstring
from hdpathmanagement import computeRelativeHTMLLink
import optiontools

#
# Module
#

def entryPoint():
    "Return information about this module to the dynamic loader."
    return {
        'name':'htmltable',
        'factory':HTMLTableFormatter,
        }



class HTMLTableFormatter(hdformatter_file.FileBasedFormatter):
    """Output documentation using HTML with tables.

      Parameters

        compactHTML -- A boolean switch to cause the formatter
                       to generate more compact HTML.  Extra
                       whitespace is removed in order to make the
                       generated files take up less space and
                       download more quickly.  The default is
                       False to cause output to be more readable.
      
        filenamePrefix -- A prefix to preprend to the base names of
                          files and directories being created.  This
                          is useful for situations where the names
                          which would be automatically generated might
                          cause a name clash or conflict.

        pageBackgroundColor -- Background color for HTML pages

        levelOneHeadingBackgroundColor -- Background color for level
                                          one heading sections.

        levelOneHeadingForegroundColor -- Foreground color for level
                                          one heading sections.

        levelTwoHeadingBackgroundColor -- Background color for level
                                          two heading sections

        levelTwoHeadingForegroundColor -- Foreground color for level
                                          two heading sections.

        dateStampFiles -- Boolean indicating whether or not to include
                          a date/time stamp in files.

        htmlQuoteText -- Boolean indicating whether or not to assume
                         that docstrings need to be quoted because
                         they might have special HTML characters in
                         them.  Defaults to true so that text is
                         quoted.

        debug -- Enable debugging comments in output.
      
    """

    def __init__(self,
                 docSet,
                 pageBackgroundColor='#ffffff',
                 levelOneHeadingBackgroundColor='#88bbee',
                 levelOneHeadingForegroundColor='#000000',
                 levelTwoHeadingBackgroundColor='#99ccff',
                 levelTwoHeadingForegroundColor='#000000',
                 appName='HappyDoc',
                 appVersion='Unknown',
                 appHome='http://happydoc.sourceforge.net',
                 styleSheet='RD.css',
                 docsetTitle=None,
                 dateStampFiles=1,
                 htmlQuoteText=1,
                 compactHTML=0,
                 debug=0,
                 **configuration):
        """Initialize the HTMLTableFormatter.

        Parameters

            'docSet' -- the DocSet instance containing global cross-reference
                      information
            
            '**configuration' -- additional, optional, configuration values

        """
        #
        # Preserve configuration parameters
        #
        self._page_bgcolor = pageBackgroundColor
        self._l1_bgcolor = levelOneHeadingBackgroundColor
        self._l1_fgcolor = levelOneHeadingForegroundColor
        self._l2_bgcolor = levelTwoHeadingBackgroundColor
        self._l2_fgcolor = levelTwoHeadingForegroundColor
        
        self._app_name = appName
        self._app_version = appVersion
        self._app_home = (appHome or 'http://www.python.org')

        self._style_sheet = styleSheet
        
        self._date_stamp_files = optiontools.getBooleanArgumentValue(dateStampFiles)
        self._html_quote_text = optiontools.getBooleanArgumentValue(htmlQuoteText)
        self._compact_html = optiontools.getBooleanArgumentValue(compactHTML)

        self.debug = debug

        #
        # Some stack counters for debugging
        #
        self._section_header_counters = {}
        self._section_header_counter = 0
        self._section_level_counter = 1
        
        #
        # Initialize the base class
        #
        apply( hdformatter_file.FileBasedFormatter.__init__,
               (self, docSet),
               configuration)
        return

    ##
    ## FileBasedFormatter implementation
    ##

    def getReference(self, infoSource, relativeSource):
        """Returns a reference to the 'infoSource' from 'relativeSource'.
        """
        #
        # Figure out the name of the infoSource
        #
        name = self.getNameForInfoSource( infoSource )
        #print '\ngetReference(', name, ',', relativeSource, ')'
        info_source_path = self.getOutputNameForObject(infoSource)
        link = computeRelativeHTMLLink(relativeSource, info_source_path,
                                       self._docset._base_directory)
        #print '>>> link to %s: %s' % (name, link)
        #if link[0] == '/':
        #    print 'STARTS AT ROOT'
            
        info = {
            'name':name,
            'link':link,
            }
        ref = '<a href="%(link)s">%(name)s</a>' % info
        return ref
    
    def getNamedReference(self, infoSource, name, relativeSource):
        """Returns a reference to 'name' within the documentation for
        'infoSource' from 'relativeSource'.
        """
        #print '\ngetNamedReference(', infoSource.getName(), ',', name, ',', relativeSource, ')'
        link = computeRelativeHTMLLink( relativeSource,
                                        self.getOutputNameForObject(infoSource),
                                        self._docset._base_directory
                                        )
        info = {
            'name':infoSource.getName(),
            'link':link,
            'target':name,
            }
        ref = '<a href="%(link)s#%(target)s">%(target)s</a>' % info
        return ref

    def getInternalReference(self, infoSource):
        """Returns a reference to 'infoSource' within the current document.
        """
        info = {
            'name':infoSource.getName(),
            }
        ref = '<a href="#%(name)s">%(name)s</a>' % info
        return ref
    
    def getPythonReference(self, moduleName):
        """Returns a reference to 'moduleName' documentation on the
        "Python.org":http://www.python.org documentation site.
        """
        if moduleName in self.sys_modules:
            return '<a href="http://www.python.org/doc/current/lib/module-%(moduleName)s.html">%(moduleName)s</a>' % locals()
        else:
            return moduleName
    
    def getFilenameExtension(self):
        "Returns the extension for creating output files."
        return 'html'

    def openOutput(self, name, title1, title2='&nbsp;'):
        """Open output destination using 'name' with the title from 'title1'.
        Write 'title2' as a secondary title to the new output.
        """
        f = hdformatter_file.FileBasedFormatter.openOutput( self,
                                                            name,
                                                            title1,
                                                            )
        self.fileHeader( title1, title2, f )
        return f

    def fileHeader(self, title1, title2='&nbsp;', output=None):
        """Write the formatting for a file header to the open file."""
        self.htmlHeader( title1, title2,
                         self._l1_bgcolor,
                         self._l1_fgcolor,
                         output) 
        return

    def closeOutput(self, output):
        "Close the 'output' handle."
        self.fileFooter(output)
        output.close()
        return

    def fileFooter(self, output):
        """Write the formatting for a file footer to the open file."""
        self.htmlFooter(output)
        return
    
    def pushSectionLevel(self, output):
        "Push a section level on the 'output' stack."
        self._section_level_counter = self._section_level_counter + 1
        self._section_header_counter = self._section_header_counters.get(
            self._section_level_counter, 0)
        self.comment('section %d:%d (push level)' % (self._section_level_counter,
                                                     self._section_header_counter),
                     output)
        self.writeHTML(
            '<table border="0" cellpadding="5" cellspacing="0" width="100%">\n',
            output)
        self.comment('push level', output)
        return

    def popSectionLevel(self, output):
        "Pop a section level from the 'output' stack."
        self.comment('section %d:%d (pop level)' % (self._section_level_counter,
                                                    self._section_header_counter),
                     output)
        #self.writeHTML('</td></tr></table>\n', output)
        self.writeHTML('</table>\n', output)
        self.comment('pop level', output)
        #
        # Depending on the pop level code to
        # close the headers for the level we just left,
        # too.
        #
        self._section_header_counters[self._section_level_counter] = 0
        #
        # Switch levels
        #
        self._section_level_counter = self._section_level_counter - 1
        #
        # Close the headers on the current level
        #
        #self._section_header_counters[self._section_level_counter] = 0
        self._section_header_counter = self._section_header_counters.get(
            self._section_level_counter, 0)
        return


    def getRootLocation(self, output,name_to_find=''):
        if name_to_find == '':
            name_to_find = self.getRootNodeName()
            
        "Return the root documentation node location relative to this 'output' location."
        first_file_opened = self.open_root_file.name
        current_output_name = output.name#[len(self._docset._base_directory)+1:]
        root_node_name = os.path.join(self._docset._base_directory,
                                      name_to_find)
        if first_file_opened == current_output_name:
            root_location = name_to_find
            #print '**SAME'
        else:
            root_location = computeRelativeHTMLLink(current_output_name,
                                                    root_node_name,
                                                    self._docset._base_directory
                                                    )
        return root_location


    def htmlHeader(self, title, subtitle, titleBg, titleFg, output):
        """Output a standard HTML header used by all output files.

        Parameters

            'title' -- title of the document

            'output' -- destination for written output

            'titleBg' -- background color for the title bar

            'titleFg' -- foreground color for text in the title bar

        """
        if not subtitle:
            subtitle = '&nbsp;'
        #
        # Determine where the root node is relative to the last
        # file opened.
        #
        root_location = self.getRootLocation(output)
        sheet_location = self.getRootLocation(output,name_to_find=self._style_sheet)
            #print '**DIFFERENT'
        #
        # Put together parts of the header
        #
        info = {
            'bgcolor':self._page_bgcolor,
            'title':title,
            'subtitle':subtitle,
            'title_bg':titleBg,
            'title_fg':titleFg,
            'root':root_location,
            'styleSheet':sheet_location,
            }
        self.writeHTML('''<html>
        
        <head>
        <title>%(title)s</title>
        <link rel="stylesheet" type="text/css" href="%(styleSheet)s">
        </head>

        <body bgcolor="%(bgcolor)s">

        <p><i><a href="%(root)s">Table of Contents</a></i></p>
        
        <table border="0" cellpadding="5" cellspacing="0" width="100%%">
        <tr bgcolor="%(title_bg)s">
            <th rowspan="2"
                valign="top"
                align="left"
                width="10%%"><font color="%(title_fg)s">%(title)s</font>
            </th>
            <th align="right"><font color="%(title_fg)s">%(subtitle)s</font></th>
        </tr>
        <tr>
        <td>
        ''' % info, output)
        self.comment('html header', output)
        return
        

    def htmlFooter(self, output):
        "Output a standard HTML footer used by all 'output' files."
        if self._date_stamp_files:
            date_str = 'on %s' % self._update_time
        else:
            date_str = ''
        info = {
            'app_name':self._app_name,
            'app_version':self._app_version,
            'app_home':self._app_home,
            'date_str':date_str,
            'root':self.getRootLocation(output),
            }
        
        self.comment('html footer', output)
        self.comment('section header %s' % str(self._section_header_counter), output)
        self.comment('section level %s' % str(self._section_level_counter), output)
        
        self.writeHTML('''
        </td>
        </tr>
        </table>

        <hr>

        <p><i><a href="%(root)s">Table of Contents</a></i></p>
        
        <p>
        This document and the software it describes are Copyright &copy; 2001-2005 Greg Landrum and Rational Discovery LLC.
        
        <p>
        <i>This document was automatically generated %(date_str)s
        by <a href="%(app_home)s">%(app_name)s</a> as modified by Greg Landrum</i>
        
        </body>
        </html>
        ''' % info, output)
        return

    def getRootNodeName(self):
        "Returns the name of the root node for documentation of this type."
        return 'index.html'

    ##
    ## HappyFormatterBase implementation
    ##

    def indent(self, output):
        "Begin an indented section."
        self.writeHTML('<ul>\n', output)
        return

    def dedent(self, output):
        "End an indented section."
        self.writeHTML('</ul>\n', output)
        return

    def writeText(self, text, output, quote=1):
        "Format and write the 'text' to the 'output'."
        if not text:
            return
        text = self._unquoteString(text)
        if self._html_quote_text and quote:
            text = StructuredText.html_quote(text)
        html = StructuredText.HTML(text)
        self.writeHTML(str(html), output)
        return

    def writeHTML(self, text, output):
        "Remove extra white space in HTML before outputting."
        #print 'text="%s"' % text
        #print 'output="%s"' % output
        compact_text = string.join( filter( None,
                                            map( string.strip,
                                                 string.split( text,
                                                               '\n'
                                                               )
                                                 )
                                            ),
                                    '\n'
                                    )
        if self._compact_html:
            self.writeRaw(compact_text, output)
        else:
            self.writeRaw(text, output)
        return

    def formatCode(self, text):
        "Format 'text' as source code and return the new string."
        formatted_text = '<pre>\n%s\n</pre>\n' % StructuredText.html_quote(text)
        return formatted_text

    def formatKeyword(self, text):
        "Format 'text' as a keyword and return the new string."
        formatted_text = '<b>%s</b>' % text
        return formatted_text

    def writeCode(self, text, output):
        "Format and write the 'text' to 'output' as source code."
        if not text:
            return
        self.writeRaw(self.formatCode(text), output)
        return

    def listHeader(self, output, title=None, allowMultiColumn=1):
        """Output 'title' as a heading for a list.  If 'allowMultiColumn' is
        true, set up the list to have more than one column.
        """
        if title:
            self.writeHTML('<h4>%s</h4>\n' % title, output)
        self.writeHTML('\n', output)
        self._pushListContext(allowMultiColumn)
        return

    def listItem(self, output, text):
        "Format and output the 'text' as a list element."
        if self.current_list_context is not None:
            self.current_list_context.append(text)
        else:
            self.writeHTML('%s<br>\n' % text, output)
        return

    def _writeListItems(self, items, output):
        "Format and output the 'items' as list elements."
        #
        # Determine the geometry of the list
        # (number of columns and rows)
        #
        num_items = len(items)
        if num_items < 10:
            num_cols = 1
            num_rows = num_items
        elif num_items < 20:
            num_cols = 2
            num_rows = 10
        else:
            num_cols = 3
            if num_items < 30:
                num_rows = 10
            else:
                num_rows = (num_items / 3) + ((num_items % 3) / 3)
        #
        # Output the list
        #
        if num_cols == 1:
            for item in self.current_list_context:
                self.writeHTML('%s<br>\n' % item, output)
        else:
            self.comment('start list', output)
            self.writeHTML('''
            <table border="0" cellspacing="2" cellpadding="2" width="100%">
              <tr>
            ''', output)
            
            for col in range(num_cols):
                self.writeHTML('<td align="LEFT" valign="TOP">\n', output)
                base_item = col * num_rows
                for item_no in range(base_item, base_item + num_rows):
                    try:
                        self.writeHTML('%s<br>\n' % items[item_no], output)
                    except IndexError:
                        break
                self.writeHTML('</td>\n', output)

            self.writeHTML('</tr>', output)
            self.writeHTML('''
            </table>
            ''', output)
            self.comment('list end', output)
        return
    
    def listFooter(self, output):
        "Write the closing footer for a list to the 'output'."
        if self.current_list_context:
            self._writeListItems(self.current_list_context, output)
        self.writeHTML('\n', output)
        self._popListContext()
        return

    def descriptiveListHeader(self, output, title):
        "Write the 'title' as the heading for a descriptive list to the 'output'."
        if title:
            self.writeHTML('<h4>%s</h4>\n' % title, output)
        self.comment('descriptive list header', output)
        self.writeHTML('<table border="0" cellpadding="5">\n', output)
        return

    def descriptiveListItem(self, output, item, description):
        "Format and write the 'item' and 'description' for a descriptive list to the 'output'."
        self.writeHTML('<tr><td valign="top" align="left">%s</td>' % item,
                      output)
        self.writeHTML('<td valign="top" align="left">', output)
        self.writeText(description, output)
        self.writeHTML('</td></tr>\n', output)
        return

    def descriptiveListFooter(self, output):
        "Write the closing footer for a descriptive list to the 'output'."
        self.writeHTML('</table>\n', output)
        self.comment('descriptive list footer', output)
        return

    def genericSectionHeader(self, output, title1, title2, anchor=None):
        """Output a standard nested table chunk which represents a section header.

        The output looks something like this::

            |--------|---------------------------|
            | title1 | title2                    |
            |        |---------------------------|
            |        | section text goes here
            |--------|

        Parameters

            'output' -- destination for written output

            'title1' -- title to be placed in left column

            'title2' -- title to be placed on top of right column

            'anchor' -- optional, anchor to which a reference can point
            to find this section

        """
        if title1 is None:
            title1 = ''
        if title2 is None:
            title2 = ''
        bgcolor = '#cccccc'
        fgcolor = '#000000'
        self._section_header_counter = self._section_header_counter + 1
        self._section_header_counters[self._section_level_counter] = self._section_header_counter
        info = {
            'bgcolor':self._l2_bgcolor,
            'fgcolor':self._l2_fgcolor,
            'title1':title1,
            'title2':title2,
            'anchor':anchor,
            }
        self.comment('section %d:%d (header)' % (self._section_level_counter,
                                                 self._section_header_counter),
                     output)
        self.writeHTML('''
        <tr>
            <th bgcolor="%(bgcolor)s"
                rowspan="2"
                valign="top"
                align="left"
                width="20%%"
                class=innerHeading
                >
                <font color="%(fgcolor)s">
                  <a name="%(anchor)s">%(title1)s</a>&nbsp;
                </font>
            </th>
            <th bgcolor="%(bgcolor)s"
                valign="top"
                align="left"
                class=innerHeading
                >
                <font color="%(fgcolor)s">%(title2)s&nbsp;</font>
            </th>
        </tr>
        <tr>
        <td>
        ''' % info, output)
        self.comment('section header', output)
        return

    def genericSectionFooter(self, output):
        "Write a general purpose section closing footer to the 'output'."
        self.comment('section %d:%d (footer)' % (self._section_level_counter,
                                                 self._section_header_counter),
                     output)
        self.writeHTML('</td></tr>\n', output)
        self.comment('section footer', output)
        self._section_header_counter = self._section_header_counter - 1
        self._section_header_counters[self._section_level_counter] = self._section_header_counter
        return

    def sectionHeader(self, output, title):
        "Write a general purpose section openning title to the 'output'."
        self.genericSectionHeader( output, title, None, title )
        return
        
    def sectionFooter(self, output):
        "Write a general purpose section closing footer to the 'output'."
        self.genericSectionFooter( output )
        return

    def itemHeader(self, output, infoObject):
        "Write a section openning header for an 'infoObject' to the 'output'."
        name = infoObject.getName()
        self.genericSectionHeader( output, None, name, name )
        return

    def itemFooter(self, output):
        "Write a section closing footer to the 'output'."
        self.genericSectionFooter( output )
        return
        
    def dividingLine(self, output, fill='-'):
        "Write a sectional dividing line made up of repeated 'fill' characters to the 'output'."
        output.write('<hr>\n')
        return

    def comment(self, text, output):
        """Output text as a comment."""
        if self.debug: self.writeHTML('<!-- %s -->\n' % text, output)
        return

    
if __name__ == '__main__':
    for fro, to in (
        ('index.html', 'index.html'),
        ('HappyDoc/CommandLineApp.py.html', 'index.html'),
        ('HappyDoc/CommandLineApp.py_CommandLineApp.html', 'index.html'),
        ('HappyDoc/CommandLineApp.py_CommandLineApp.html',
         'HappyDoc/CommandLineApp.py.html'),
        ('/home/hellmann/devel/HappyDoc/doc/HappyDoc/ts_regex.py_compile.html',
         '/home/hellmann/devel/HappyDoc/doc/index.html'),
        
        ):
        #print 'FROM: ', fro
        #print 'TO  : ', to
        #print 'LINK: ', computeRelativeHTMLLink(fro, to)
        computeRelativeHTMLLink(fro, to)
        print
