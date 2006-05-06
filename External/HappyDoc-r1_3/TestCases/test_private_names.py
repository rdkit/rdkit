#!/usr/bin/env python
#
# $Id: test_private_names.py,v 1.1 2001/03/21 21:00:20 glandrum Exp $
#
# Time-stamp: <01/02/03 12:56:57 dhellmann>
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

"""Module with some public and some private names.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: test_private_names.py,v $',
    'rcs_id'       : '$Id: test_private_names.py,v 1.1 2001/03/21 21:00:20 glandrum Exp $',
    'creator'      : 'Doug Hellmann <DougHellmann@bigfoot.com>',
    'project'      : 'HappyDoc',
    'created'      : 'Sat, 03-Feb-2001 12:56:49 EST',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 21:00:20 $',
}

#
# Import system modules
#


#
# Import Local modules
#


#
# Module
#

class Public:
    """Test case for public/private names.
    
    If HappyDoc works correctly, you should only see the
    definition of 'public_method' below when the '--no_private_names'
    option is specified.
    """

    def public_method(self):
        "This method is public."
        pass

    def _private_method(self):
        "This method is private."
        pass

    def __getattr__(self, attrname):
        "This method should not be hidden."
        pass

    
