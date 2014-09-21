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
from xml.etree import ElementTree
# check the version of ElementTree.  We need at least version 1.2
# in order for the XPath-style parsing stuff to work
import re
vers = re.split("[a-zA-Z]",ElementTree.VERSION)[0]
if vers < '1.2':
  raise ImportError('The PubMed record interface requires a version of ElementTree >= 1.2')


class Record(object):
  def __init__(self,element):
    for field in self._fieldsOfInterest:
        setattr(self,field,'')
    self._element = element
  def toXML(self):
    from cStringIO import StringIO
    sio = StringIO()
    ElementTree.ElementTree(self._element).write(sio)
    return sio.getvalue()
  
class SummaryRecord(Record):
  _fieldsOfInterest=['PubMedId','PubDate','Source','Authors',
                    'Title','Volume','Issue','Pages','Lang',
                    'HasAbstract','RecordStatus']
  def __init__(self,element):
    Record.__init__(self,element)
    for item in element.getiterator('Item'):
      if item.attrib['Name'] in self._fieldsOfInterest:
        setattr(self,item.attrib['Name'],item.text)
    if self.PubDate:
      self.PubYear = str(self.PubDate).split(' ')[0]

class JournalArticleRecord(Record):
  _fieldsOfInterest=['PubMedId','PubYear','Source','Authors',
                    'Title','Volume','Issue','Pages','Lang',
                    'Abstract']
  def __init__(self,element):
    Record.__init__(self,element)

    cite = self._element.find('MedlineCitation')
    self.PubMedId = cite.findtext('PMID')
    article = cite.find('Article')
    issue = article.find('Journal/JournalIssue')
    self.Volume = issue.findtext('Volume')
    self.Issue = issue.findtext('Issue')
    self.PubYear = issue.findtext('PubDate/Year')
    if not self.PubYear:
      txt = issue.findtext('PubDate/MedlineDate')
      self.PubYear = txt.split(' ')[0]
    self.Title = unicode(article.findtext('ArticleTitle'))
    self.Pages = article.findtext('Pagination/MedlinePgn')
    abs = article.findtext('Abstract/AbstractText')
    if abs:
      self.Abstract = unicode(abs)

    self.authors = []
    tmp = []
    for author in article.find('AuthorList').getiterator('Author'):
      last = unicode(author.findtext('LastName'))
      first = unicode(author.findtext('ForeName'))
      initials =  unicode(author.findtext('Initials'))
      self.authors.append((last,first,initials))
      tmp.append('%s %s'%(last,initials))
    self.Authors=', '.join(tmp)
    journal = cite.findtext('MedlineJournalInfo/MedlineTA')
    if journal:
      self.Source = unicode(journal)

    self.ParseKeywords()
    self.ParseChemicals()

  def ParseKeywords(self):
    self.keywords = []
    headings = self.find('MedlineCitation/MeshHeadingList')
    if headings:
      for heading in headings.getiterator('MeshHeading'):
        kw = unicode(heading.findtext('DescriptorName'))
        for qualifier in heading.getiterator('QualifierName'):
          kw += ' / %s'%(unicode(qualifier.text))
        self.keywords.append(kw)
    
  def ParseChemicals(self):
    self.chemicals = []
    chemicals = self.find('MedlineCitation/ChemicalList')
    if chemicals:
      for chemical in chemicals.getiterator('Chemical'):
        name = chemical.findtext('NameOfSubstance').encode('utf-8')
        rn = chemical.findtext('RegistryNumber').encode('utf-8')
        if rn != '0':
          self.chemicals.append('%s <%s>'%(name,rn))
        else:
          self.chemicals.append('%s'%(name))
    

  # --------------------------------------------
  #
  #  We'll expose these ElementTree methods in case
  #  client code wants to pull extra info
  #
  def getiterator(self,key=None):
    if key is not None:
      return self._element.getiterator(key)
    else:
      return self._element.getiterator()
  def find(self,key):
    return self._element.find(key)
  def findtext(self,key):
    return self._element.findtext(key)
  def findall(self,key):
    return self._element.findall(key)
  
class LinkRecord(Record):
  _fieldsOfInterest=[]
  def __init__(self,element):
    Record.__init__(self,element)
    self.PubMedId = self._element.text
    nbr = self._element.get('HasNeighbor','N')
    if nbr == 'Y':
      self.HasNeighbor = 1
    else:
      self.HasNeighbor = 0
      
    
    
