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
import copy

urlBase="http://eutils.ncbi.nlm.nih.gov/entrez/eutils"
searchBase=urlBase+"/esearch.fcgi"
summaryBase=urlBase+"/esummary.fcgi"
fetchBase=urlBase+"/efetch.fcgi"
linkBase=urlBase+"/elink.fcgi"
# for links to pubmed web pages:
queryBase="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi"

_details = {
  'db':'pubmed',
  'retmode':'xml',
  'tool':'RDPMTools',
  'email':'Info@RationalDiscovery.com',
  }

def details():
  return copy.copy(_details)

# FIX: Allow PMID searches
searchableFields={
  "Author":("AU","Authors' Name "),
  "Keyword":("MH","MeSH keyword"),
  "Substance":("NM","Substance Name"),
  "Title":("TI","Text from the article title"),
  "Title/Abstract":("TIAB","Text from the article title and/or abstract"),
  "Registry Number":("RN","CAS Registry Number"),
  "Subject":("SB","Pubmed Subject Subset (tox,aids,cancer,bioethics,cam,history,space,systematic)"),
  "Journal":("TA","Journal Name"),
  "Year":("DP","Publication Date"),
  "Affiliation":("AD","Authors' affiliation"),
  }
searchableFieldsOrder=[
  "Author",
  "Keyword",
  "Title",
  "Title/Abstract",
  "Substance",
  "Registry Number",
  "Subject",
  "Journal",
  "Year",
  "Affiliation",
  ]
