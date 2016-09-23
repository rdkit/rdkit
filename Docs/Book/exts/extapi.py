# Original code from: http://git.savannah.gnu.org/cgit/kenozooid.git/tree/doc/extapi.py
#
# Kenozooid - software stack to support different capabilities of dive
# computers.
#
# Copyright (C) 2009 by Artur Wroblewski <wrobell@pld-linux.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import os.path
from docutils import nodes


def api_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
  """
    Role `:api:` bridges generated API documentation by tool like EpyDoc
    with Sphinx Python Documentation Generator.

    Other tools, other than EpyDoc, can be easily supported as well.

    First generate the documentation to be referenced, i.e. with EpyDoc::

        $ mkdir -p doc/_build/html/api
        $ epydoc -o doc/_build/html/api ...

    Next step is to generate documentation with Sphinx::

        $ sphinx-build doc doc/_build/html

    """
  basedir = 'api'
  prefix = '_build/html/'  # fixme: fetch it from configuration
  exists = lambda f: os.path.exists(prefix + f)

  # assume module is references
  name = '%s' % text
  uri = file = '%s/%s-module.html' % (basedir, text)
  chunks = text.split('.')

  # if not module, then a class
  if not exists(file):
    name = text.split('.')[-1]
    uri = file = '%s/%s-class.html' % (basedir, text)

  # if not a class, then function or class method
  if not exists(file):
    method = chunks[-1]
    fprefix = '.'.join(chunks[:-1])
    # assume function is referenced
    file = '%s/%s-module.html' % (basedir, fprefix)
    if exists(file):
      uri = '%s#%s' % (file, method)
    else:
      # class method is references
      file = '%s/%s-class.html' % (basedir, fprefix)
      if exists(file):
        name = '.'.join(chunks[-2:])  # name should be Class.method
        uri = '%s/%s-class.html#%s' % (basedir, fprefix, method)

  if exists(file):
    node = nodes.reference(rawtext, name, refuri=uri, **options)
  else:
    # cannot find reference, then just inline the text
    node = nodes.literal(rawtext, text)
  return [node], []


def setup(app):
  app.add_role('api', api_role)
