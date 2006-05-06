//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

/*! \file RDKitBase.h

  \brief pulls in the core \c RDKit functionality

*/

#ifndef _RDKIT_BASE_H
#define _RDKIT_BASE_H

#include <RDGeneral/Invariant.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SanitException.h>

#if 0
namespace RDKit{
  template class AtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class AtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class AromaticAtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class AromaticAtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class HeteroatomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class HeteroatomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
  template class QueryAtomIterator_<Atom,ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::type>;
  template class QueryAtomIterator_<const Atom,const ROMol,ROMol::GRAPH_MOL_ATOM_PMAP::const_type>;
};
#endif


#endif
