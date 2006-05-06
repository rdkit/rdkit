#!/usr/bin/env python
#
# Time-stamp: <01/02/03 12:52:34 dhellmann>
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


"""Documentation is directed to stdout.

$Id: docset_stdout.py,v 1.1 2001/03/21 21:00:00 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: docset_stdout.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 18:27:39 EDT',
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
import sys

#
# Import Local modules
#
import docset_singlefile

#
# Module
#

def entryPoint():
    "Return info about this module to the dynamic loader."
    return { 'name':'stdout',
             'factory':StdOutDocSet,
             }



class StdOutDocSet(docset_singlefile.SingleFileDocSet):
    """Documentation is directed to stdout.

      Parameters

        *See DocSet.*

    """

    def openOutput(self, name, title, subtitle=''):
        "Returns sys.stdout as the output destination, always."
        stdout = sys.stdout
        stdout.write('%s\n' % title)
        if subtitle:
            stdout.write('  %s\n' % subtitle)
        stdout.write('\n')
        return stdout

    def closeOutput(self, output):
        """No-op"""
        return

    def close(self):
        """No-op"""
        return
