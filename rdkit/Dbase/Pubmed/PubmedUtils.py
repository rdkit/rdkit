# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from __future__ import print_function
import re,urllib
import QueryParams


def FillParagraph(text,width=70,indentLinesBy=0):
  """

  >>> txt = "123 567 90 123\\n456 89\\n"
  >>> print FillParagraph(txt,width=10)
  123 567 90
  123 456 89
  
  >>> print FillParagraph(txt,width=10,indentLinesBy=3)
  123 567 90
     123 456
     89
  >>> txt = '123 567 \\n'
  >>> print FillParagraph(txt,width=10)
  123 567

  """
  expr = re.compile('[\ \t\n\r]+')
  words = expr.split(text)
  res = []
  line = ''
  for word in words:
    if word:
      if len(line)+len(word)>=width:
        res.append(line)
        line = ' '*indentLinesBy
        line += word
      else:
        if line:
          line += ' %s'%word
        else:
          line = word
  if line:
    res.append(line)    
  return '\n'.join(res)

def RecordsToPubmedText(records):
  """  Exports a sequence of JournalArticleRecords to a Medline text
  format (not XML) that can be imported by EndNote

  """
  res = []
  for i in range(len(records)):
    rec = records[i]
    res.append('UI  - %d'%i)
    res.append('PMID- %s'%rec.PubMedId)
    for auth in rec.authors:
      res.append('AU  - %s %s'%(auth[0],auth[2]))
    res.append('TI  - %s'%FillParagraph(rec.Title,indentLinesBy=6))
    res.append('TA  - %s'%rec.Source)
    res.append('DP  - %s'%rec.PubYear)
    res.append('PG  - %s'%rec.Pages)
    res.append('VI  - %s'%rec.Volume)
    res.append('IP  - %s'%rec.Issue)
    res.append('AB  - %s'%FillParagraph(rec.Abstract,indentLinesBy=6))

    for keyword in rec.keywords:
      res.append('MH  - %s'%(keyword))
                 
  return '\n'.join(res)
    


def GetSearchURL(term,db="PubMed"):
  queryParams = {'CMD':'search',
                 'DB':db,
                 'term':term}
  url = '%s?%s'%(QueryParams.queryBase,
                 urllib.urlencode(queryParams))
  return url
  


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys,os.path
  failed,tried = _test()
  sys.exit(failed)

  import Searches
  recs = Searches.GetRecords(['11960484','10893315'],
                             conn=open('test_data/records.xml','r'))
  res = RecordsToPubmedText(recs)
  print(res)
