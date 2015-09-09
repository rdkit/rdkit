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

%template(FilterCatalogEntry_Vect) std::vector< boost::shared_ptr<RDKit::FilterCatalogEntry> >;

%template(FilterCatalogEntryVect) std::vector< const RDKit::FilterCatalogEntry* >;

%include "enums.swg"

%include <../RDGeneral/Dict.h>
%include <../Catalogs/Catalog.h>
%include <../Catalogs/CatalogParams.h>
%include <GraphMol/Substruct/SubstructMatch.h>


%newobject RDKit::FilterCatalogEntry::getProp;
%extend RDKit::FilterCatalogEntry {
  std::string getProp(const std::string key){
    std::string res;
    ($self)->getProp(key, res);
    return res;
  }
 }



%extend RDKit::FilterCatalog {
  bool canSerialize() const {
    return RDKit::FilterCatalogCanSerialize();
  }

  // boost does a bad job of wrapping shared_ptr<const T> so we will do the
  //  unthinkable and cast away const.
  //  Also we can't use FilterCatalog::SENTRY because swig thinks it is a new
  //   type.  Bad swig!
  boost::shared_ptr<RDKit::FilterCatalogEntry> getFirstMatch(const ROMol &mol) {
    RDKit::FilterCatalog::CONST_SENTRY res = self->getFirstMatch(mol);
    return boost::const_pointer_cast<RDKit::FilterCatalogEntry>(res);
  }
  
  std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry> > getMatches(const ROMol &mol) {
    std::vector<RDKit::FilterCatalog::CONST_SENTRY> matches = self->getMatches(mol);
    std::vector<RDKit::FilterCatalog::SENTRY> res;
    res.reserve(matches.size());
    for (size_t i=0; i< matches.size(); ++i) {
      res.push_back( boost::const_pointer_cast<RDKit::FilterCatalogEntry>(matches[i]) );
    }
    return res;
 }
}

%ignore RDKit::FilterCatalog::getFirstMatch;
%ignore RDKit::FilterCatalog::getMatches;
//%ignore RDKit::FilterCatalogEntry::getPropList;
%ignore RDKit::Dict::getPropList;

%include <GraphMol/FilterCatalog/FilterMatcherBase.h>
%include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
%include <GraphMol/FilterCatalog/FilterCatalog.h>



