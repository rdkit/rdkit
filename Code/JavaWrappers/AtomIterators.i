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

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"

%{
#include <GraphMol/AtomIterators.h>
#include <GraphMol/ROMol.h>
%}


%include <GraphMol/AtomIterators.h>
/* Equality testing operators have been overloaded, so we need to expose them in a different way.
   These 'rename' lines appear to have to precede the 'template' definitions. */
%rename(eq) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>::operator==;
%rename(ne) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>::operator!=;
%rename(eq) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>::operator==;
%rename(ne) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>::operator!=;
%rename(eq) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator==;
%rename(ne) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator!=;
%rename(eq) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator==;
%rename(ne) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator!=;
/* Increment and decrement operators currently necessary */
%rename(next) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>::operator++;
%rename(next) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>::operator++;
%rename(next) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator++;
%rename(next) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator++;
%rename(prev) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>::operator--;
%rename(prev) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>::operator--;
%rename(prev) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator--;
%rename(prev) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator--;
/* A better name for the iterator's Atom object */
%rename(getAtom) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>::operator*;
%rename(getAtom) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>::operator*;
%rename(getAtom) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator*;
%rename(getAtom) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>::operator*;

%template(AtomIterator) RDKit::AtomIterator_<RDKit::Atom,RDKit::ROMol>;
%template(HeteroatomIterator) RDKit::HeteroatomIterator_<RDKit::Atom,RDKit::ROMol>;
%template(AromaticAtomIterator) RDKit::AromaticAtomIterator_<RDKit::Atom,RDKit::ROMol>;
%template(QueryAtomIterator) RDKit::QueryAtomIterator_<RDKit::Atom,RDKit::ROMol>;



