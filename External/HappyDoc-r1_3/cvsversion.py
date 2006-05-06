#!/usr/bin/env python
#
# $Id: cvsversion.py,v 1.1 2001/03/21 20:59:26 glandrum Exp $
#
# Time-stamp: <01/02/03 12:49:12 dhellmann>
#
# Copyright 2001 Doug Hellmann.
#
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

"""Get the CVS version information based on the $Name:  $ token.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: cvsversion.py,v $',
    'rcs_id'       : '$Id: cvsversion.py,v 1.1 2001/03/21 20:59:26 glandrum Exp $',
    'creator'      : 'Doug Hellmann <DougHellmann@bigfoot.com>',
    'project'      : 'HappyDoc',
    'created'      : 'Sat, 03-Feb-2001 12:48:26 EST',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 20:59:26 $',
}

#
# Import system modules
#
import string

#
# Import Local modules
#


#
# Module
#

def cvs_product_version():
    cvs_version_string='$Name:  $'
    cvs_version_parts=string.split(cvs_version_string)
    if len(cvs_version_parts) >= 3:
        app_version = string.strip(cvs_version_parts[1])
    else:
        app_version = 'WORKING'
    return app_version

if __name__ == '__main__':
    print cvs_product_version()
    
