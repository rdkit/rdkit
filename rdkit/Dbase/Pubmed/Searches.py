# $Id$
#
# Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Tools for doing PubMed searches and processing the results

NOTE: much of the example code in the documentation here uses XML
files from the test_data directory in order to avoid having to call
out to PubMed itself.  Actual calls to the functions would not include
the _conn_ argument.

"""
from __future__ import print_function
from rdkit import RDConfig
import QueryParams,Records
import urllib,urllib2
from xml.etree import ElementTree

def openURL(url,args):
  proxy_support = urllib2.ProxyHandler({})
  opener = urllib2.build_opener(proxy_support)
  conn = urllib2.urlopen(url,args)
  return conn
  
def GetNumHits(query,url=QueryParams.searchBase):
  """ returns a tuple of pubmed ids (strings) for the query provided

  To do a search, we need a query object:
  >>> query = QueryParams.details()

  set up the search parameters:
  >>> query['term'] = 'penzotti je AND grootenhuis pd'
  >>> query['field'] = 'auth'

  now get the search ids:
  >>> counts = GetNumHits(query)
  >>> counts
  2

  alternately, we can search using field specifiers:
  >>> query = QueryParams.details()
  >>> query['term'] = 'penzotti je[au] AND hydrogen bonding[mh]'
  >>> counts = GetNumHits(query)
  >>> counts
  3
  

  """
  query['rettype']='count'
  conn = openURL(url,urllib.urlencode(query))
  pubmed = ElementTree.parse(conn)
  countText = pubmed.findtext('Count')
  if countText:
    res = int(countText)
  else:
    res = 0
  return res


def GetSearchIds(query,url=QueryParams.searchBase):
  """ returns a tuple of pubmed ids (strings) for the query provided

  To do a search, we need a query object:
  >>> query = QueryParams.details()

  set up the search parameters:
  >>> query['term'] = 'penzotti je AND grootenhuis pd'
  >>> query['field'] = 'auth'

  now get the search ids:
  >>> ids = GetSearchIds(query)
  >>> len(ids)
  2
  >>> ids[0]
  '11960484'
  >>> ids[1]
  '10893315'
  

  """
  conn = openURL(url,urllib.urlencode(query))
  pubmed = ElementTree.parse(conn)
  res = [id.text for id in pubmed.getiterator('Id')]
  return tuple(res)

def GetSummaries(ids,query=None,url=QueryParams.summaryBase,conn=None):
  """ gets a set of document summary records for the ids provided

  >>> ids = ['11960484']
  >>> summs = GetSummaries(ids,conn=open(os.path.join(testDataDir,'summary.xml'),'r'))
  >>> len(summs)
  1
  >>> rec = summs[0]
  >>> isinstance(rec,Records.SummaryRecord)
  1
  >>> rec.PubMedId
  '11960484'
  >>> rec.Authors
  'Penzotti JE, Lamb ML, Evensen E, Grootenhuis PD'
  >>> rec.Title
  'A computational ensemble pharmacophore model for identifying substrates of P-glycoprotein.'
  >>> rec.Source
  'J Med Chem'
  >>> rec.Volume
  '45'
  >>> rec.Pages
  '1737-40'
  >>> rec.HasAbstract
  '1'

  """
  if not conn:
    try:
      iter(ids)
    except TypeError:
      ids = [ids,]
    if not query:
      query = QueryParams.details()
    ids = map(str,ids)
    query['id'] = ','.join(ids)
    conn = openURL(url,urllib.urlencode(query))
  pubmed = ElementTree.parse(conn)
  res = []
  for summary in pubmed.getiterator('DocSum'):
    rec = Records.SummaryRecord(summary)
    if rec.PubMedId in ids:
      res.append(rec)
      ids.remove(rec.PubMedId)

  return tuple(res)

def GetRecords(ids,query=None,url=QueryParams.fetchBase,conn=None):
  """ gets a set of document summary records for the ids provided

  >>> ids = ['11960484']
  >>> recs = GetRecords(ids,conn=open(os.path.join(testDataDir,'records.xml'),'r'))
  >>> len(recs)
  1
  >>> rec = recs[0]
  >>> rec.PubMedId
  '11960484'
  >>> rec.Authors
  u'Penzotti JE, Lamb ML, Evensen E, Grootenhuis PD'
  >>> rec.Title
  u'A computational ensemble pharmacophore model for identifying substrates of P-glycoprotein.'
  >>> rec.Source
  u'J Med Chem'
  >>> rec.Volume
  '45'
  >>> rec.Pages
  '1737-40'
  >>> rec.PubYear
  '2002'
  >>> rec.Abstract[:10]
  u'P-glycopro'

  We've also got access to keywords:
  >>> str(rec.keywords[0])
  'Combinatorial Chemistry Techniques'
  >>> str(rec.keywords[3])
  'Indinavir / chemistry'

  and chemicals:
  >>> rec.chemicals[0]
  'P-Glycoprotein'
  >>> rec.chemicals[2]
  'Nicardipine <55985-32-5>'
  
  
  """
  if not conn:
    try:
      iter(ids)
    except TypeError:
      ids = [ids,]
    if not query:
      query = QueryParams.details()
    query['id'] = ','.join(map(str,ids))
    conn = openURL(url,urllib.urlencode(query))

  pubmed = ElementTree.parse(conn)
  res = []
  for article in pubmed.getiterator('PubmedArticle'):
    rec = Records.JournalArticleRecord(article)
    if rec.PubMedId in ids:
      res.append(rec)
  return tuple(res)

def CheckForLinks(ids,query=None,url=QueryParams.linkBase,conn=None):
  if not conn:
    try:
      iter(ids)
    except TypeError:
      ids = [ids,]
    if not query:
      query = QueryParams.details()
    query['id'] = ','.join(map(str,ids))
    conn = openURL(url,urllib.urlencode(query))
    query['cmd'] = 'ncheck'
  pubmed = ElementTree.parse(conn)

  checklist = pubmed.find('LinkSet/IdCheckList')
  recs = [Records.LinkRecord(x) for x in checklist.getiterator('Id')]

  res = {}
  for rec in recs:
    id = rec.PubMedId
    res[id] = rec.HasNeighbor
  return res

def GetLinks(ids,query=None,url=QueryParams.linkBase,conn=None):
  if not conn:
    try:
      iter(ids)
    except TypeError:
      ids = [ids,]
    if not query:
      query = QueryParams.details()
    query['id'] = ','.join(map(str,ids))
    conn = openURL(url,urllib.urlencode(query))
    query['cmd'] = 'neighbor'

  pubmed = ElementTree.parse(conn)
  linkset = pubmed.find('LinkSet/LinkSetDb')
  scores = []
  scoreNorm = 1.0
  for link in linkset.getiterator('Link'):
    id = link.findtext('Id')
    score = float(link.findtext('Score'))
    scores.append([id,score])
    # we'll normalize scores by the score for the first of the query ids:
    if id == ids[0]:
      scoreNorm = score
  for i in range(len(scores)):
    id,score = scores[i]
    scores[i] = id,score/scoreNorm
  return tuple(scores)


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys,os.path
  testDataDir = os.path.join(RDConfig.RDCodeDir,'Dbase','Pubmed','test_data')
  failed,tried = _test()
  sys.exit(failed)
  #query = QueryParams.details()
  #query['term']='landrum ga'
  #query['field']='auth'
  #ids = GetSearchIds(query)
  #print ids
  #ids = ids[:2]
  ids = ['11666868','11169640']
  if 0:
    summs = GetSummaries(ids,conn=open('summary.xml','r'))
    print('summs:',summs)
    for summary in summs:
      print(summary.Authors)
      print('\t',summary.Title)
      print('\t',summary.Source,end='')
      print(summary.Volume,end='')
      print(summary.Pages,end='')
      print(summary.PubDate)
    
  if 0:
    ids = ['11666868']
    res = GetRecords(ids,conn=open('records.xml','r'))
    for record in res:
      print(record.Authors)
      print('\t',record.Title)
      print('\t',record.Journal,end='')
      print(record.Volume,end='')
      print(record.Pages,end='')
      print(record.PubYear)
      print()
    
  if 0:
    ids = ['11666868','11169640']
    res = CheckForLinks(ids,conn=open('haslinks.xml','r'))
    print(res)
    
  if 0:
    ids = ['11666868']
    res = GetLinks(ids,conn=open('links.xml','r'))
    #res = GetLinks(ids)
    for id,score in res[:10]:
      print(id,score)

