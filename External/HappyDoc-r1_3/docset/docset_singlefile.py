#!/usr/bin/env python
#
# Time-stamp: <01/02/03 12:52:46 dhellmann>
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


"""All of the documentation goes into one file.

$Id: docset_singlefile.py,v 1.1 2001/03/21 21:00:00 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: docset_singlefile.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 18:27:00 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 21:00:00 $',
    }

#
# Import system modules
#


#
# Import Local modules
#
import docset_multiplefile

#
# Module
#

def entryPoint():
    "Return info about this module to the dynamic loader."
    return { 'name':'singlefile',
             'factory':SingleFileDocSet,
             }



class SingleFileDocSet(docset_multiplefile.DocSet):
    """All of the documentation is in one file.

      Parameters

        *See DocSet.*
        
    """

    def _getOutputHandle(self, name, title, subtitle):
        return docset_multiplefile.DocSet.openOutput( self,
                                                      name,
                                                      title,
                                                      subtitle)        
    
    def openOutput(self, name, title, subtitle=''):
        """Return a handle to the single output file."""
        #
        # If we already have an open handle, return it
        #
        if hasattr(self, '_output_handle'):
            self._formatter.comment('OPEN: %s' % name, self._output_handle)
            self._formatter.sectionHeader(self._output_handle, title)
            self._output_handle_open_count = self._output_handle_open_count + 1
            return self._output_handle
        #
        # This assumes that we are called first with a reasonable name
        # for the single output file.
        #
        output_handle = self._getOutputHandle(name, title, subtitle)
        #
        # Store the handle for use later.
        #
        self._output_handle = output_handle
        self._output_handle_open_count = 1
        return output_handle

    def closeOutput(self, output):
        """No-op"""
        self._output_handle_open_count = self._output_handle_open_count - 1
        self._formatter.sectionFooter(output)
        return

    def close(self):
        self._formatter.popSectionLevel(self._output_handle)
        docset_multiplefile.DocSet.closeOutput(self, self._output_handle)
        return
    
    def _writeTOC(self):
        """No-op"""
        return
    
    def _describeClassInModuleNode(self, output, class_output_file_name, class_info):
        return
    
