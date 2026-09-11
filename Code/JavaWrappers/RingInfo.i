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
#include <map>
#include <vector>
#include <GraphMol/RingInfo.h>
%}

// Keep the managed-language API based on owned vectors. The C++ methods now
// return non-owning views whose lifetimes cannot be represented safely by the
// Java and C# wrappers.
%ignore RDKit::RingInfo::RingsView;
%ignore RDKit::RingInfo::atomRings;
%ignore RDKit::RingInfo::bondRings;
%ignore RDKit::RingInfo::atomMembers;
%ignore RDKit::RingInfo::bondMembers;
%ignore RDKit::RingInfo::addRing;
%rename(atomRings) RDKit::RingInfo::atomRingsForWrapper;
%rename(bondRings) RDKit::RingInfo::bondRingsForWrapper;
%rename(atomMembers) RDKit::RingInfo::atomMembersForWrapper;
%rename(bondMembers) RDKit::RingInfo::bondMembersForWrapper;
%rename(addRing) RDKit::RingInfo::addRingForWrapper;

%include <GraphMol/RingInfo.h>

%extend RDKit::RingInfo {
  RDKit::VECT_INT_VECT atomRingsForWrapper() const {
    RDKit::VECT_INT_VECT result;
    const auto rings = $self->atomRings();
    result.reserve(rings.size());
    for (const auto ring : rings) {
      result.emplace_back(ring.begin(), ring.end());
    }
    return result;
  }

  RDKit::VECT_INT_VECT bondRingsForWrapper() const {
    RDKit::VECT_INT_VECT result;
    const auto rings = $self->bondRings();
    result.reserve(rings.size());
    for (const auto ring : rings) {
      result.emplace_back(ring.begin(), ring.end());
    }
    return result;
  }

  RDKit::INT_VECT atomMembersForWrapper(unsigned int idx) const {
    const auto members = $self->atomMembers(idx);
    return {members.begin(), members.end()};
  }

  RDKit::INT_VECT bondMembersForWrapper(unsigned int idx) const {
    const auto members = $self->bondMembers(idx);
    return {members.begin(), members.end()};
  }

  unsigned int addRingForWrapper(const RDKit::INT_VECT &atomIndices,
                                 const RDKit::INT_VECT &bondIndices) {
    return $self->addRing(atomIndices, bondIndices);
  }
}
