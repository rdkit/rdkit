#
#  Copyright (C) 2004  Rational Discovery LLC
#   All Rights Reserved
#
""" utility functions for lists

"""


def CompactListRepr(lst):
  """

  >>> CompactListRepr([0,1,1,1,1,0])
  '[0]+[1]*4+[0]'
  >>> CompactListRepr([0,1,1,2,1,1])
  '[0]+[1]*2+[2]+[1]*2'
  >>> CompactListRepr([])
  '[]'
  >>> CompactListRepr((0,1,1,1,1))
  '[0]+[1]*4'
  >>> CompactListRepr('foo')
  "['f']+['o']*2"

  """
  if not len(lst):
    return '[]'

  components = []
  last = lst[0]
  count = 1
  i = 1
  while i < len(lst):
    if lst[i] != last:
      label = '[%s]' % repr(last)
      if count > 1:
        label += '*%d' % count
      components.append(label)
      count = 1
      last = lst[i]
    else:
      count += 1
    i += 1
  if count != 0:
    label = '[%s]' % repr(last)
    if count > 1:
      label += '*%d' % count
    components.append(label)

  return '+'.join(components)


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
