/* 
 *  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
 *  All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: 
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following 
 *       disclaimer in the documentation and/or other materials provided 
 *       with the distribution.
 *     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
 *       nor the names of its contributors may be used to endorse or promote 
 *       products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//%import "ROMol.i"
%include "std_string.i"
%include "std_vector.i"
%{
#include <../RDGeneral/Dict.h>
#include <../Catalogs/Catalog.h>
#include <../Catalogs/CatalogParams.h>
#include <GraphMol/FilterCatalog/FilterMatcherBase.h>
#include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
  
  // bug fix for swig, it removes these from their namespaces
  typedef RDCatalog::Catalog<RDKit::FilterCatalogEntry, RDKit::FilterCatalogParams>::paramType_t paramType_t;
  typedef RDCatalog::Catalog<RDKit::FilterCatalogEntry, RDKit::FilterCatalogParams>::entryType_t entryType_t;
  typedef std::vector<std::string> STR_VECT;

%}

%include "enums.swg"
%include <boost_shared_ptr.i>

%include <../RDGeneral/Dict.h>
%include <../Catalogs/Catalog.h>
%include <../Catalogs/CatalogParams.h>
%include <GraphMol/Substruct/SubstructMatch.h>

%template(SENTRY) boost::shared_ptr<RDKit::FilterCatalogEntry>;
%template(FilterCatalogEntryVect) std::vector< const RDKit::FilterCatalogEntry* >;

%newobject RDKit::FilterCatalogEntry::getProp;
%extend RDKit::FilterCatalogEntry {
  std::string getProp(const std::string key){
    std::string res;
    ($self)->getProp(key, res);
    return res;
  }
 }



%extend RDKit::FilterCatalog {
  // not a new object, don't delete from the outside
  const RDKit::FilterCatalogEntry * getFirstMatch(const ROMol &mol) const {
    // dangerous, this can be delete from underneath you
    return self->getFirstMatch(mol).get();
  }

  std::vector<const RDKit::FilterCatalogEntry *> getMatches(const ROMol &mol) const {
    // dangerous, this can be delete from underneath you

    std::vector<boost::shared_ptr<const RDKit::FilterCatalogEntry> > results = \
      self->getMatches(mol);
    
    std::vector<const RDKit::FilterCatalogEntry*> ptrs;
    ptrs.reserve(results.size());
    for(size_t i=0;i<results.size();++i) {
      ptrs.push_back(results[i].get());
    }
    return ptrs;
  }
  
  bool canSerialize() const {
    return RDKit::FilterCatalogCanSerialize();
  }
 }

%ignore RDKit::FilterCatalog::getFirstMatch;
%ignore RDKit::FilterCatalog::getMatches;
//%ignore RDKit::FilterCatalogEntry::getPropList;
%ignore RDKit::Dict::getPropList;

%include <GraphMol/FilterCatalog/FilterMatcherBase.h>
%include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
%include <GraphMol/FilterCatalog/FilterCatalog.h>



