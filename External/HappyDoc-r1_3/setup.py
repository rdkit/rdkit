#!/usr/bin/env python
#
# $Id: setup.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $
#
# Time-stamp: <01/02/03 12:51:39 dhellmann>
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

"""Distutils setup file for HappyDoc

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: setup.py,v $',
    'rcs_id'       : '$Id: setup.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $',
    'creator'      : 'Doug Hellmann <DougHellmann@bigfoot.com>',
    'project'      : 'UNSPECIFIED',
    'created'      : 'Sat, 03-Feb-2001 12:51:26 EST',

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
from distutils.core import setup

#
# Import Local modules
#
from cvsversion import cvs_product_version

#
# Module
#

setup ( name = "HappyDoc",
        version = cvs_product_version(),
        description = "Python HappyDoc ",
        author = "Doug Hellmann",
        author_email = "doughellmann@bigfoot.com",
        url = "http://sourceforge.net/projects/happydoc",
        packages = [ '',
                     'docset',
                     'hdformatter',
                     ],
        package_dir = { '': '.' },
        extra_path= 'happydoc',
        scripts = ['happydoc'],
        )

