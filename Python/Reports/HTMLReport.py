#
#  Copyright (C) 2003-2005  Rational Discovery LLC
#   All Rights Reserved
#
import RDConfig

def FormatTable(data,rowHeadings=0,colHeadings=0,border=0,extraParams=''):
  """ data should be a sequence of sequences

  >>> d = ((1,2),(3,4))
  >>> print FormatTable(d)
  <table border=0 >
    <tr>
      <td>1</td>
      <td>2</td>
    </tr>
    <tr>
      <td>3</td>
      <td>4</td>
    </tr>
  </table>
  >>> print FormatTable(d,rowHeadings=1,border=1)
  <table border=1 >
    <tr>
      <th>1</th>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
    </tr>
  </table>
  >>> print FormatTable(d,colHeadings=1)
  <table border=0 >
    <tr>
      <th>1</th>
      <th>2</th>
    </tr>
    <tr>
      <td>3</td>
      <td>4</td>
    </tr>
  </table>
  

  """
  if not (data and hasattr(data,'__len__') and len(data)):
    return ''

  lines = ['<table border=%(border)d %(extraParams)s>'%locals()]
  nRows = len(data)
  nCols = len(data[0])
  for i in range(nRows):
    row = data[i]
    lines.append('  <tr>')
    for j in range(nCols):
      cell = row[j]
      if colHeadings and i==0:
        lines.append('    <th>%s</th>'%str(cell))
      elif rowHeadings and j==0:
        lines.append('    <th>%s</th>'%str(cell))
      else:
        lines.append('    <td>%s</td>'%str(cell))
    lines.append('  </tr>')
  lines.append('</table>')
  res = '\n'.join(lines)
  return res
  

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)

