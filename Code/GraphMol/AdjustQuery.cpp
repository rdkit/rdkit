//
//  Copyright (c) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

#include <vector>
#include <algorithm>

namespace RDKit {
namespace {
bool isMapped(const Atom *atom) {
  return atom->hasProp(common_properties::molAtomMapNumber);
}
}

namespace MolOps {
ROMol *adjustQueryProperties(const ROMol &mol,
                             const AdjustQueryParameters *params) {
  RWMol *res = new RWMol(mol);
  try {
    adjustQueryProperties(*res, params);
  } catch (MolSanitizeException &se) {
    delete res;
    throw se;
  }
  return static_cast<ROMol *>(res);
}
void adjustQueryProperties(RWMol &mol, const AdjustQueryParameters *inParams) {
  AdjustQueryParameters params;
  if (inParams) {
    params = *inParams;
  }
  const RingInfo *ringInfo = mol.getRingInfo();

  if (params.aromatizeIfPossible) {
    unsigned int failed;
    sanitizeMol(mol, failed, SANITIZE_SYMMRINGS | SANITIZE_SETAROMATICITY);
  } else {
    if (!ringInfo->isInitialized()) {
      MolOps::symmetrizeSSSR(mol);
    }
  }

  if (params.makeAtomsGeneric) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!((params.makeAtomsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numAtomRings(i)) &&
          !((params.makeAtomsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numAtomRings(i)) &&
          !((params.adjustDegreeFlags & ADJUST_IGNOREMAPPED) &&
            isMapped(mol.getAtomWithIdx(i)))) {
        QueryAtom *qa = new QueryAtom();
        qa->setQuery(makeAtomNullQuery());
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
      }
    }
  }  // end of makeAtomsGeneric
  if (params.makeBondsGeneric) {
    for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
      if (!((params.makeBondsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numBondRings(i)) &&
          !((params.makeBondsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numBondRings(i))) {
        QueryBond *qb = new QueryBond();
        qb->setQuery(makeBondNullQuery());
        const bool preserveProps = true;        
        mol.replaceBond(i, qb, preserveProps);
        delete qb;
      }
    }
  }  // end of makeBondsGeneric
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    Atom *at = mol.getAtomWithIdx(i);
    // pull properties we need from the atom here, once we
    // create a query atom they may no longer be valid.
    unsigned int nRings = ringInfo->numAtomRings(i);
    int atomicNum = at->getAtomicNum();
    if (params.makeDummiesQueries && atomicNum == 0 && !at->hasQuery() &&
        !at->getIsotope()) {
      QueryAtom *qa = new QueryAtom();
      qa->setQuery(makeAtomNullQuery());
      const bool updateLabel = false;
      const bool preserveProps = true;
      mol.replaceAtom(i, qa, updateLabel, preserveProps);
      delete qa;
      at = mol.getAtomWithIdx(i);
    }  // end of makeDummiesQueries
    if (params.adjustDegree &&
        !((params.adjustDegreeFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustDegreeFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustDegreeFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) &&
        !((params.adjustDegreeFlags & ADJUST_IGNORENONDUMMIES) && atomicNum) &&
        !((params.adjustDegreeFlags & ADJUST_IGNOREMAPPED) && isMapped(at))) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        qa = new QueryAtom(*at);
        const bool updateLabel = false;        
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      if (params.adjustDegree == 2) {
        ATOM_EQUALS_QUERY *tmp = new ATOM_EQUALS_QUERY;
        tmp->setVal(qa->getTotalDegree() - qa->getTotalNumHs(true));
        tmp->setDataFunc(queryAtomHeavyAtomDegree);
        tmp->setDescription("queryAtomHeavyDegree");
        qa->expandQuery(tmp);
      } else {
        qa->expandQuery(makeAtomExplicitDegreeQuery(qa->getDegree()));
      }
    }  // end of adjust degree
    if (params.adjustRingCount &&
        !((params.adjustRingCountFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustRingCountFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustRingCountFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) &&
        !((params.adjustRingCountFlags & ADJUST_IGNORENONDUMMIES) &&
          !((params.adjustRingCountFlags & ADJUST_IGNOREMAPPED) &&
            isMapped(at)) &&
          atomicNum)) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        qa = new QueryAtom(*at);
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      qa->expandQuery(makeAtomInNRingsQuery(nRings));
    }  // end of adjust ring count
  }    // end of loop over atoms
  if (params.makeBondsGeneric) {
    ROMol::EDGE_ITER firstB, lastB;
    boost::tie(firstB, lastB) = mol.getEdges();
    while (firstB != lastB) {
      BOND_SPTR bond = mol[*firstB];
      ++firstB;
    }
  }
}
}  // end of MolOps namespace
}  // end of RDKit namespace
