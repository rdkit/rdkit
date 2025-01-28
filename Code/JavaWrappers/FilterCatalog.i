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
%include <boost_shared_ptr.i>
%shared_ptr(RDKit::FilterCatalogEntry)

%{
#include <../RDGeneral/Dict.h>
#include <../Catalogs/Catalog.h>
#include <../Catalogs/CatalogParams.h>
#include <GraphMol/FilterCatalog/FilterMatcherBase.h>
#include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>

// bug fix for swig, it removes these from their namespaces
typedef RDCatalog::Catalog<RDKit::FilterCatalogEntry, RDKit::FilterCatalogParams>::paramType_t paramType_t;

typedef std::vector<std::string> STR_VECT;
%}

%template(FilterCatalogEntry_Vect) std::vector< boost::shared_ptr<RDKit::FilterCatalogEntry> >;
%template(FilterCatalogEntry_VectVect) std::vector<std::vector< boost::shared_ptr<const RDKit::FilterCatalogEntry> > >;

%template(FilterCatalogEntryVect) std::vector< const RDKit::FilterCatalogEntry* >;
%template(UChar_Vect) std::vector<unsigned char>;
%template(FilterMatch_Vect) std::vector<RDKit::FilterMatch>;
%include "enums.swg"

%include <../RDGeneral/Dict.h>
%include <../Catalogs/Catalog.h>
%include <../Catalogs/CatalogParams.h>
%include <GraphMol/Substruct/SubstructMatch.h>


%typemap(javacode) RDKit::FilterCatalog %{
  public static FilterCatalog Deserialize(byte[] b) {
    UChar_Vect vec = null;
    try {
      vec = new UChar_Vect();
      vec.reserve(b.length);
      for (int size=0;size<b.length;++size) {
        vec.add((short)b[size]);
      }
      return new FilterCatalog(vec);
    } finally {
      if (vec != null) {
        vec.delete();
      }
    }
  }
%}

%extend RDKit::FilterMatch {
  MatchVectType getAtomMatches() {
    return ($self)->atomPairs;
  }
 }

%newobject RDKit::FilterCatalogEntry::getProp;
%extend RDKit::FilterCatalogEntry {
  std::string getProp(const std::string key){
    std::string res;
    ($self)->getProp(key, res);
    return res;
  }

  std::vector<RDKit::FilterMatch> getFilterMatches(const ROMol &mol) {
    std::vector<RDKit::FilterMatch> matches;
    ($self)->getFilterMatches(mol, matches);
    return matches;
  }
 }

%extend RDKit::FilterCatalog {
  FilterCatalog(const std::vector<unsigned char> & data ) {
    std::string str(data.begin(), data.end());
    return new RDKit::FilterCatalog(str);
  }

  static bool canSerialize() {
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

  // re-wrap swig is making duplicate entries for some strange reason
  unsigned int addEntry(boost::shared_ptr<RDKit::FilterCatalogEntry> entry) {
    return self->addEntry(entry);
  }

  bool removeEntry(boost::shared_ptr<RDKit::FilterCatalogEntry> entry) {
    return self->removeEntry(entry);
  }

  // swig const-ptr shenanigans again
  boost::shared_ptr<RDKit::FilterCatalogEntry> getEntry(unsigned int idx) const {
    return boost::const_pointer_cast<RDKit::FilterCatalogEntry>(
      self->getEntry(idx));
  }

  unsigned int getIdxForEntry(const boost::shared_ptr<FilterCatalogEntry> &entry) const {
    return self->getIdxForEntry(entry);
  }
}

%ignore RDKit::FilterCatalog(const std::string &);
%ignore RDKit::FilterCatalog::getFirstMatch;
%ignore RDKit::FilterCatalog::getMatches;
%ignore RDKit::FilterCatalog::addEntry;
%ignore RDKit::FilterCatalog::removeEntry;
%ignore RDKit::FilterCatalog::setCatalogParams;
%ignore RDKit::FilterCatalog::getIdxForEntry;
%ignore RDKit::FilterCatalog::getEntryWithIdx;

%ignore RDKit::FilterMatch::operator==;
%ignore std::vector<RDKit::FilterMatch>::operator==;

//%ignore RDKit::FilterCatalogEntry::getPropList;
%ignore RDKit::Dict::getPropList;


#ifdef SWIGJAVA
%typemap(jni) std::string RDKit::FilterCatalog::Serialize "jbyteArray"
%typemap(jtype) std::string RDKit::FilterCatalog::Serialize "byte[]"
%typemap(jstype) std::string RDKit::FilterCatalog::Serialize "byte[]"
%typemap(javaout) std::string RDKit::FilterCatalog::Serialize {
  return $jnicall;
}
%typemap(out) std::string RDKit::FilterCatalog::Serialize {
  $result = JCALL1(NewByteArray, jenv, $1.size());
  JCALL4(SetByteArrayRegion, jenv, $result, 0, $1.size(), (const jbyte*)$1.c_str());
}
#endif

#ifdef SWIGCSHARP
%typemap(csbase) RDKit::FilterCatalogParams::FilterCatalogs "uint"
#endif

%immutable RDKit::DEFAULT_FILTERMATCHERBASE_NAME;

%include <GraphMol/FilterCatalog/FilterMatcherBase.h>
%include <GraphMol/FilterCatalog/FilterCatalogEntry.h>
%include <GraphMol/FilterCatalog/FilterCatalog.h>



%pragma(java) modulecode=%{
   public static FilterCatalog FilterCatalogDeserialize(byte[] b) {
     UChar_Vect vec = new UChar_Vect();
     vec.reserve(b.length);
     for (int size=0;size<b.length;++size) {
       vec.add((short)b[size]);
     }
     return new FilterCatalog(vec);
   }
%}
