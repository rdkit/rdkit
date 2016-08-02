//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include "Utilites.h"
#include "Pattern.h"

namespace RDKit {
 namespace StructureCheck {

static inline
RDKit::Bond::BondType convertBondType(AABondType bt) {
    static const RDKit::Bond::BondType rdbt[] = {
        RDKit::Bond::UNSPECIFIED,  //BT_NONE = 0,
        RDKit::Bond::SINGLE,
        RDKit::Bond::DOUBLE,
        RDKit::Bond::TRIPLE,
        RDKit::Bond::AROMATIC,
        RDKit::Bond::UNSPECIFIED,   // ??? SINGLE_DOUBLE = 5,
        RDKit::Bond::ONEANDAHALF,   //SINGLE_AROMATIC = 6,
        RDKit::Bond::TWOANDAHALF,   //DOUBLE_AROMATIC = 7,
    };
    return bt<= DOUBLE_AROMATIC ? rdbt[bt] : RDKit::Bond::OTHER;    // ??
}

static inline
AABondType convertBondType(RDKit::Bond::BondType rdbt) {
    static const AABondType bt[] = {
        BT_NONE,
        SINGLE,
        DOUBLE,
        TRIPLE,
        BT_NONE,    // QUADRUPLE,
        BT_NONE,    // QUINTUPLE,
        BT_NONE,    // HEXTUPLE,
        SINGLE_AROMATIC,    // ONEANDAHALF
        DOUBLE_AROMATIC,    // TWOANDAHALF
        BT_NONE,    // THREEANDAHALF,
        BT_NONE,    // FOURANDAHALF,
        BT_NONE,    // FIVEANDAHALF,
        AROMATIC,
    };
    return rdbt <= RDKit::Bond::AROMATIC ? bt[rdbt] : BT_NONE; //??
}


/*
* Checks if the ligand is compatible with the bond and the atom.
*/
bool LigandMatches(const Atom &a, const Bond &b, const Ligand &l, bool use_charge) {
    if (l.BondType != ANY_BOND  && convertBondType(b.getBondType()) != l.BondType)
        return false;

    if (l.Radical != ANY_RADICAL && a.getNumRadicalElectrons() != l.Radical)
        return false;
    if ((l.Charge != ANY_CHARGE || use_charge)
        && l.Charge != ANY_CHARGE  && a.getFormalCharge != l.Charge)
        return false;
    return (AtomSymbolMatch(a.getSymbol(), l.AtomSymbol));
}

struct AtomNeighbor {
    unsigned AtomIdx;   // the molecule's atom
    unsigned BondIdx;   // the molecule's bond
    AtomNeighbor (unsigned a, unsigned b) : AtomIdx(a), BondIdx(b) {}
    };

bool TransformAugmentedAtoms(RWMol &mol, const std::vector<std::pair<AugmentedAtom, AugmentedAtom> > &aapair) {
    bool transformed = false;
// find neighbors:
    std::vector< std::vector<AtomNeighbor> > neighbors(mol.getNumAtoms());
    std::vector<bool> bondInRing(mol.getNumBonds());
    std::vector<bool> atomInRing(mol.getNumAtoms());
    std::vector<std::pair<unsigned, unsigned> > bondsToRemove;

    for (unsigned i = 0; i < mol.getNumBonds(); i++) {
        const Bond* bond = mol.getBondWithIdx(i);
        unsigned a1 = bond->getEndAtomIdx();
        unsigned a2 = bond->getBeginAtomIdx();
        neighbors[a1].push_back(AtomNeighbor(a2, i));
        neighbors[a2].push_back(AtomNeighbor(a1, i));
        bondInRing[i] = false;                      //init
        atomInRing[a1] = atomInRing[a2] = false;    //init
    }

    const RingInfo::VECT_INT_VECT& rings = mol.getRingInfo()->bondRings();
    for (RingInfo::VECT_INT_VECT::const_iterator r = rings.begin(); r != rings.end(); r++) {
        for (size_t j = 0; j < r->size(); j++) {
            bondInRing[(*r)[j]] = true;
            const Bond* bond = mol.getBondWithIdx((*r)[j]);
            unsigned a1 = bond->getEndAtomIdx();
            unsigned a2 = bond->getBeginAtomIdx();
            atomInRing[a1] = atomInRing[a2] = true;
        }
    }
// find matched atoms and replace
    for (size_t i = 0; i < aapair.size(); i++) {
        //AAFix():
        for (unsigned j = 0; j < mol.getNumAtoms(); j++) {
            const Atom* atom = mol.getAtomWithIdx(j);
            const AugmentedAtom &aa1 = aapair[i].first;
            // check if atom match to aapair[i].first atom and topology
            if (neighbors[j].size() == aa1.Ligands.size()
                && (ANY_CHARGE == aa1.Charge  || atom->getFormalCharge() == aa1.Charge)
                && (ANY_RADICAL== aa1.Radical || atom->getNumRadicalElectrons() == aa1.Radical)
                && ( TP_NONE == aa1.Topology
                  ||(RING  == aa1.Topology && atomInRing[j])
                  ||(CHAIN == aa1.Topology && ! atomInRing[j]))
                && AtomSymbolMatch(atom->getSymbol(), aa1.AtomSymbol)) {
                 //RecMatch():
                unsigned matched = 0;
                const AugmentedAtom &aa2 = aapair[i].second;
                std::vector<bool> visited(aa2.Ligands.size());
                std::vector<unsigned> match(aa2.Ligands.size());
                for (size_t l = 0; l < aa2.Ligands.size(); l++) {
                    visited[l] = false;
                    match  [l] = -1;
                }

                for (size_t k = 0; k < neighbors[j].size(); k++) {
                    const Atom* nbrAtom = mol.getAtomWithIdx(neighbors[j][k].AtomIdx);
                    const Bond* nbrBond = mol.getBondWithIdx(neighbors[j][k].BondIdx);
                    for (size_t l = 0; l < aa2.Ligands.size(); l++)
                     if ( ! visited[l]) {
                        const Ligand &ligand = aa2.Ligands[l];
                        if ((0 == ligand.SubstitutionCount || neighbors[j].size() == ligand.SubstitutionCount)
                            && (ANY_CHARGE == ligand.Charge  || nbrAtom->getFormalCharge() == ligand.Charge)
                            && (ANY_RADICAL== ligand.Radical || nbrAtom->getNumRadicalElectrons() == ligand.Radical)
                            &&((ANY_BOND  == ligand.BondType || convertBondType(nbrBond->getBondType()) == ligand.BondType)
                             ||(SINGLE_DOUBLE == ligand.BondType 
                                && (RDKit::Bond::SINGLE == nbrBond->getBondType() 
                                 || RDKit::Bond::DOUBLE == nbrBond->getBondType())))
                            && AtomSymbolMatch(nbrAtom->getSymbol(), ligand.AtomSymbol)) {
                            matched++;
                            match  [l] = neighbors[j][k].AtomIdx;
                            visited[l] = true;
                            break;
                        }
                    }
                }
                if(matched != neighbors[j].size())
                    continue; // go to next molecule's atom

                // Replace Augmented Atom. TransformAA()
                // change central atom
                Atom* a = mol.getAtomWithIdx(j);
                if (a->getSymbol() != aa2.AtomSymbol)
                    a->setAtomicNum(getAtomicNumber(aa2.AtomSymbol));
                if (a->getFormalCharge() != aa2.Charge)
                    a->setFormalCharge(aa2.Charge);
                if (a->getNumRadicalElectrons() != aa2.Radical)
                    a->setNumRadicalElectrons(aa2.Radical);
                // change ligand atoms and their bonds with central atom
                for (size_t l = 0; l < aa2.Ligands.size(); l++) {
                    Atom* al = mol.getAtomWithIdx(match[l]);
                    const Ligand &ligand = aa2.Ligands[l];
                    if (al->getSymbol() != ligand.AtomSymbol)
                        al->setAtomicNum(getAtomicNumber(ligand.AtomSymbol));
                    if (al->getFormalCharge() != ligand.Charge)
                        al->setFormalCharge(ligand.Charge);
                    if (al->getNumRadicalElectrons() != ligand.Radical)
                        al->setNumRadicalElectrons(ligand.Radical);

                    if (BT_NONE == aa2.Ligands[l].BondType) // remove bond
                        bondsToRemove.push_back(std::pair<unsigned, unsigned>(j, neighbors[j][l].AtomIdx));
                    else if (aa1.Ligands[l].BondType != aa2.Ligands[l].BondType)
                            mol.getBondWithIdx(neighbors[j][l].BondIdx)->setBondType(
                                convertBondType(aa2.Ligands[l].BondType));
                }
                transformed = true;
            }
        }
    }
    // remove bonds marked to delete // BT_NONE
    if (!bondsToRemove.empty()) {
        for (size_t i = 0; i < bondsToRemove.size(); i++)
            mol.removeBond(bondsToRemove[i].first, bondsToRemove[i].second);
    }
    return transformed;
}

//-----------------------------------------------------------------------------
// CheckAtoms():
//-----------------------------------------------------------------------------

/*
* Recursively searches *mp for the next atom to be appended to
* match[0..level]. Returns TRUE if all ligands in aap are mapped
* successfully. nbp[i] describes the neighbour bonds and atoms
* of atom i+1.
*/
static
bool RecMatch(const ROMol &mol, std::vector<unsigned> &match, unsigned int level,
              const AugmentedAtom &aa,  const std::vector<Neighbourhood> &nbp) {
    bool is_new;

    if (level == aa.Ligands.size())
        return true;
    const Neighbourhood &nbph = nbp[match[0]];
    for (unsigned i = 0; i < nbph.Atoms.size(); i++) {
        const Atom &atom = *mol.getAtomWithIdx(nbph.Atoms[i]);
        const Bond &bond = *mol.getBondWithIdx(nbph.Bonds[i]);
        if ((atom.getFormalCharge() == aa.Ligands[level].Charge || aa.Ligands[level].Charge == ANY_CHARGE)
            && (atom.getNumRadicalElectrons() == aa.Ligands[level].Radical || aa.Ligands[level].Radical == ANY_RADICAL)
            && convertBondType(bond.getBondType()) == aa.Ligands[level].BondType
            && (aa.Ligands[level].SubstitutionCount == 0 ||
                nbp[nbph.Atoms[i]].Atoms.size() == aa.Ligands[level].SubstitutionCount)
            && AtomSymbolMatch(atom.getSymbol(), aa.Ligands[level].AtomSymbol)) {
            is_new = true;
            match[level + 1] = nbph.Atoms[i];
            for (unsigned j = 0; j <= level; j++)
                if (nbph.Atoms[i] == match[j])
                    is_new = false;
            if (is_new && RecMatch(mol, match, level + 1, aa, nbp))
                return true;
        }
    }
    return false;
}

/*
* Tests if atom i in *mp matches the augmented atom description
* *aap. nbp[] is used to speed up access to neighbour atoms and
* bonds. The first matching atom mapping is placed into match[1..].
* i is stored in match[0].
*/
static
bool AAMatch(const ROMol &mol, unsigned i,
    std::vector<unsigned> &match,
    const AugmentedAtom &aa,
    const std::vector<unsigned> &atom_ring_status,
    const std::vector<Neighbourhood> &nbp) {

    const Atom &atom = *mol.getAtomWithIdx(i);
    if (nbp[i].Atoms.size() == aa.Ligands.size() 
     && (aa.Charge == ANY_CHARGE || atom.getFormalCharge() == aa.Charge)
     && (aa.Radical == ANY_RADICAL || atom.getNumRadicalElectrons() == aa.Radical) &&
        AtomSymbolMatch(atom.getSymbol(), aa.AtomSymbol))
    {
        if ( ! atom_ring_status.empty() && aa.Topology == RING   &&  atom_ring_status[i] == 0)
            return false;
        if ( ! atom_ring_status.empty() && aa.Topology == CHAIN  &&  atom_ring_status[i] != 0)
            return false;
        match[0] = i;
        return RecMatch(mol, match, 0, aa, nbp);
    }
    else 
        return false;
}


/*
* Computes how many basis rings each bond shares and how many
* ring bonds are attached to an atom. The results are stored in
* atom_status[] and bond_status[] respectively.
*/
static
void RingState(const ROMol & mol, std::vector<unsigned>& atom_status, std::vector<unsigned> &bond_status) {
        atom_status.resize(mol.getNumAtoms());
        bond_status.resize(mol.getNumBonds());
        for (unsigned i = 0; i < atom_status.size(); i++)
            atom_status[i] = 0;
        for (unsigned i = 0; i < bond_status.size(); i++)
            bond_status[i] = 0;
        if (mol.getNumBonds() == 0)
            return;
        // for each bond compute amount of rings that contains the bond
        const RDKit::RingInfo &ringInfo = *mol.getRingInfo();
        const VECT_INT_VECT &bondRings = ringInfo.bondRings();
        for (unsigned ir = 0; ir < bondRings.size(); ir++) {
            const INT_VECT &ring = bondRings[ir];
            for (unsigned i = 0; i < ring.size(); i++) {
                bond_status[i]++;
            }
        }
        for (unsigned i = 0; i < bond_status.size(); i++)
            if (bond_status[i] > 0) {
                const Bond &bond = *mol.getBondWithIdx(i);
                atom_status[bond.getBeginAtomIdx()] ++;
                atom_status[bond.getEndAtomIdx()  ] ++;
            }
}

/*
* Checks if every atom in *mp matches one of the augmented atoms
* in good_atoms[0..ngood-1]. It returns TRUE if all atoms gave a match or FALSE otherwise.
*/
bool CheckAtoms(const ROMol &mol, const std::vector<AugmentedAtom> &good_atoms) {
    std::vector<Neighbourhood> neighbours(mol.getNumAtoms());
    std::vector<unsigned> match;    //[MAXNEIGHBOURS + 1];
    std::vector<unsigned> atom_status(mol.getNumAtoms());
    std::vector<unsigned> bond_status(mol.getNumBonds());

    RingState(mol, atom_status, bond_status);
    SetupNeighbourhood(mol, neighbours);

    unsigned nmatch = 0;
    if ( ! good_atoms.empty())
     for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
        for (unsigned j = 0; j < good_atoms.size(); j++)
        {
            // check for ring state of central atom. Ring matches to ring only
            if (good_atoms[j].Topology == RING   && 0 == atom_status[i])
                continue;
            if (good_atoms[j].Topology == CHAIN  && 0 != atom_status[i])
                continue;
            if (neighbours[i].Atoms.size() == good_atoms[j].Ligands.size()
             && AAMatch(mol, i, match, good_atoms[j], atom_status, neighbours)) {
                nmatch++;
                break;
            }
        }
    }
    return (nmatch == mol.getNumAtoms());
}

 }// namespace StructureCheck
} // namespace RDKit
