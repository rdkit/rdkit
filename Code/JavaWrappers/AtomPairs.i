/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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

%{
#include <boost/cstdint.hpp>
#include <DataStructs/BitVects.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
%}
%include <DataStructs/BitVects.h>
%include <DataStructs/SparseIntVect.h>

%rename(eq) RDKit::SparseIntVect<boost::uint32_t>::operator==;
%rename(ne) RDKit::SparseIntVect<boost::uint32_t>::operator!=;
%rename(eq) RDKit::SparseIntVect<boost::int32_t>::operator==;
%rename(ne) RDKit::SparseIntVect<boost::int32_t>::operator!=;
%rename(eq) RDKit::SparseIntVect<boost::int64_t>::operator==;
%rename(ne) RDKit::SparseIntVect<boost::int64_t>::operator!=;

%template(SparseIntVectu32) RDKit::SparseIntVect<boost::uint32_t>;
%template(SparseIntVect32) RDKit::SparseIntVect<boost::int32_t>;
%template(SparseIntVect64) RDKit::SparseIntVect<boost::int64_t>;

%ignore RDKit::SparseIntVect<boost::uint32_t>::getNonzeroElements const;
%ignore RDKit::SparseIntVect<boost::int32_t>::getNonzeroElements const;
%ignore RDKit::SparseIntVect<boost::int64_t>::getNonzeroElements const;

%extend RDKit::SparseIntVect<boost::uint32_t> {
  std::vector<std::pair<boost::uint32_t, int> > getNonzero() const{
    std::vector<std::pair<boost::uint32_t, int> > res;
    for(std::map<boost::uint32_t,int>::const_iterator es=$self->getNonzeroElements().begin();
        es!=$self->getNonzeroElements().end();++es){
      res.push_back(std::make_pair((boost::uint32_t)es->first,(int)es->second));
    }
    return res;
  }
}
%extend RDKit::SparseIntVect<boost::int32_t> {
  std::vector<std::pair<boost::int32_t, int> > getNonzero() const{
    std::vector<std::pair<boost::int32_t, int> > res;
    for(std::map<boost::int32_t,int>::const_iterator es=$self->getNonzeroElements().begin();
        es!=$self->getNonzeroElements().end();++es){
      res.push_back(std::make_pair((boost::int32_t)es->first,(int)es->second));
    }
    return res;
  }
}
%extend RDKit::SparseIntVect<boost::int64_t> {
  std::vector<std::pair<boost::int64_t, int> > getNonzero() const{
    std::vector<std::pair<boost::int64_t, int> > res;
    for(std::map<boost::int64_t,int>::const_iterator es=$self->getNonzeroElements().begin();
        es!=$self->getNonzeroElements().end();++es){
      res.push_back(std::make_pair((boost::int64_t)es->first,(int)es->second));
    }
    return res;
  }
}
%newobject getAtomPairFingerprint;
%newobject getHashedAtomPairFingerprint;
%newobject getHashedAtomPairFingerprintAsBitVect;
%newobject getTopologicalTorsionFingerprint;
%newobject getHashedTopologicalTorsionFingerprint;
%newobject getHashedTopologicalTorsionFingerprintAsBitVect;
%include <GraphMol/Fingerprints/AtomPairs.h>
