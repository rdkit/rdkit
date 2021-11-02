//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Stereo.h"
#include "Pattern.h"
#include "ReCharge.h"

namespace RDKit {
namespace StructureCheck {

/*
 * Returns the total charge of all atoms in molecule.
 */
int TotalCharge(const ROMol &mol) {
  int charge = 0;
  for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
    const Atom &atom = *mol.getAtomWithIdx(i);
    charge += atom.getFormalCharge();
  }
  return 0;
}

/*
 * Checks if all defined stereocenters in *mp are unambiguously
 * numbered by numbering, i.e. all four ligands get different
 * numbers. It returns TRUE if this is the case and FALSE otherwise.
 */
bool AllCentersRefined(const ROMol &mol, std::vector<unsigned> &numbering) {
  std::vector<Neighbourhood> neighbour_array;
  int parity;
  unsigned refnum[4], nref;

  SetupNeighbourhood(mol, neighbour_array);
  for (unsigned i = 0; i < neighbour_array.size(); i++) {
    unsigned element = mol.getAtomWithIdx(i)->getAtomicNum();
    unsigned n_ligands = (unsigned)neighbour_array[i].Bonds.size();
    if ((6 == element       // "C"
         || 16 == element   // "S"
         || 7 == element    // "O"
         || 8 == element    // "N"
         || 15 == element   // "P"
         || 14 == element)  // "Si"
        && n_ligands > 2 && n_ligands <= 4) {
      const Neighbourhood &nbp = neighbour_array[i];
      parity = AtomParity(mol, i, nbp);  // avalon: i + 1 ???
      if (parity == EVEN_PARITY || parity == ODD_PARITY) {
        // extraxt numbers of stereoligands
        nref = (unsigned)nbp.Atoms.size();
        for (unsigned j = 0; j < nref; j++) refnum[j] = numbering[nbp.Atoms[j]];
        // sort ligands
        for (unsigned j = nref; j > 0; j--)
          for (unsigned k = j; k > 0; k--)
            if (refnum[k] < refnum[k - 1]) {
              unsigned h = refnum[k];
              refnum[k] = refnum[k - 1];
              refnum[k - 1] = h;
            }
        // check if there is a duplicate
        for (unsigned j = 1; j < nref; j++)
          if (refnum[j] == refnum[j - 1]) {
            // sprintf(msg_buffer, "%10s    unrefined parity at atom %d",
            // mp->name, i + 1);
            // AddMsgToList(msg_buffer);
            return false;
          }
      }
    }
  }
  return true;
}

/*
 * Removes hydrogens from *mp until desired_charge is reached. The
 * positions for hydrogen removal are selected by "acidity" combined
 * with a refinement algorithm. It returns TRUE if molecule could be
 * neutralized and FALSE if any problem were encountered.
 * *ndeprot and *nrefine are set to the number of deprotonations
 * and refinement cycles performed.
 */
bool ChargeFix::rechargeMolecule(unsigned &ndeprot, unsigned &nrefine) {
  bool changed = false;
  int nacid = 0;
  double gap = 0.0, pKa_value = 0.0;
  double acidity_limit = 24.0;
  std::vector<unsigned> numbering(Mol.getNumAtoms());

  ndeprot = 0;  // number of deprotonation cycles
  nrefine = 0;  // number of refinements necessary
  for (;;) {
    while (TotalCharge(Mol) > Options.DesiredCharge) {
      setpKaValues();
      nacid = markMostAcidicAtoms(gap, acidity_limit);
      if (pKa_value > acidity_limit) {
        //                    sprintf(msg_buffer, "pKa_value = %.2g",
        //                    pKa_value);
        //                    AddMsgToList(msg_buffer);
        return false;  // not acidic enough
      } else if (nacid == 1 ||
                 (nacid == TotalCharge(Mol) && gap > 8.0)) {  // easy case
        // sprintf(msg_buffer, "pKa = %.2g", pKa_value);
        // AddMsgToList(msg_buffer);
        decrementMarkedCharges();
      } else
        break;
      ndeprot++;
    }

    if (TotalCharge(Mol) > Options.DesiredCharge) {
      nrefine++;
      //// ?? expl. H ??                numbering = TypeAlloc(mp->n_atoms, int);
      nacid = refineAcidicAtoms(numbering);
      if (nacid == 1)
        decrementMarkedCharges();
      else if (AllCentersRefined(Mol, numbering)) {
        for (unsigned i = 0; i < Mol.getNumAtoms(); i++)  // Select one mark
          if (AtomColor[i] != 0 && nacid != 1) {
            nacid--;
            AtomColor[i] = 0;
          }
        decrementMarkedCharges();
      } else {
        //                    sprintf(msg_buffer, "%10s: could not fix charges",
        //                    mp->name);
        //                    AddMsgToList(msg_buffer);
        return false;
      }
      ndeprot++;
    } else {
      resetValues();
      return true;
    }
  }
  return changed;
}

void ChargeFix::resetColors() {
  BondColor.resize(Mol.getNumBonds());
  AtomColor.resize(Mol.getNumAtoms());
  for (unsigned int &i : BondColor) i = 0;
  for (unsigned int &i : AtomColor) i = 0;
}

void ChargeFix::resetValues() {
  AtompKaValue.resize(Mol.getNumAtoms());
  for (unsigned i = 0; i < AtomColor.size(); i++) AtompKaValue[i] = 0.0;
}

/*
 * Estimates the pKa value of acidic atoms in molecule *mp and sets
 * the 'value' fields to this value.
 * It returns TRUE if all went well and FALSE otherwise.
 */
bool ChargeFix::setpKaValues() {
  const ROMol &mol = Mol;
  bool result = true;
  const unsigned na = mol.getNumAtoms();
  std::vector<Neighbourhood> neighbour_array(na);
  std::vector<unsigned> match;  //[MAXNEIGHBOURS + 1];

  // const AugmentedAtom *AAp;

  //    struct inc_entry_t  *cip, *ip;
  //    struct path_entry_t *ppa, *ppb;

  //    struct reaccs_atom_t *aap, *bap;
  //    struct reaccs_bond_t *abp, *bbp;

  SetupNeighbourhood(Mol, neighbour_array);

  resetColors();
  AtomOldpKaValue = AtompKaValue;
  resetValues();

  for (unsigned i = 0; i < na; i++) {
    const Atom *ap = Mol.getAtomWithIdx(i);
    // Is atom acidic?
    bool found = false;
    const std::vector<unsigned> dummy_atom_ring_status;
    for (const auto &AcidicAtom : Options.AcidicAtoms)
      if (AAMatch(Mol, i, AcidicAtom, dummy_atom_ring_status, neighbour_array,
                  Options.Verbose)) {
        AtomColor[i]++;
        AtompKaValue[i] = 0.0;
        // AAp = &AcidicAtom;
        found = true;
        break;
      }
    if (!found) continue;  // center not acidic -> goto next atom

    //        if (charge_log && old_values[i] != 0.0) // Prints the header of an
    //        atom prediction line.
    //            StartPredictionLine(AAp->ShortName, i + 1, old_values[i]);
    //        if (old_values[i] != 0.0)
    //            fprintf(stderr, "atom %d: '%s'\n", i + 1, AAp->short_name); //
    //            throw..

    // add local charge increment
    for (const auto &j : Options.ChargeIncTable)
      if (AtomSymbolMatch(ap->getSymbol(), j.AtomSymbol)) {
        AtompKaValue[i] += ap->getFormalCharge() * j.LocalInc;
        //                if (charge_log && old_values[i] != 0 && ap->charge !=
        //                NONE) Options.ChargeIncTable[j].local_inc_used++;
        //                if (charge_log && old_values[i] != 0)
        //                fprintf(charge_log, "+%d*%c3", ap->charge, 'B' + j);
        break;
      }

    // add local atom acidity (atom_acidity_table)
    for (const auto &j : Options.AtomAcidity)
      if (AtomSymbolMatch(ap->getSymbol(), j.AtomSymbol)) {
        AtompKaValue[i] += j.LocalInc;  // ip->LocalInc
        //                if (charge_log && old_values[i] != 0)
        //                ip->local_inc_used++;
        //                if (charge_log && old_values[i] != 0)
        //                    fprintf(charge_log, "+%c10", 'B' + j);
        break;
      }

    const Neighbourhood &nbp = neighbour_array[i];  // for all alpha neighbours
    for (unsigned j = 0; j < nbp.Atoms.size(); j++) {
      const Atom &aap = *Mol.getAtomWithIdx(nbp.Atoms[j]);  // alpha atom
      const Bond &abp = *Mol.getBondWithIdx(nbp.Bonds[j]);  // alpha bond
      const PathEntry *ppa = nullptr;
      // fetch alpha Path conductivity
      found = false;
      for (const auto &k : Options.AlphaPathTable)
        if (AtomSymbolMatch(ap->getSymbol(), k.Path.AtomSymbol) &&
            LigandMatches(aap, abp, k.Path.Ligands[0], false)) {
          ppa = &k;
          found = true;
          break;
        }
      if (!found) {
        // sprintf(msg_buffer, "%10s: no alpha Path for atom %d", mp->name, i +
        // 1);
        // AddMsgToList(msg_buffer);
        result = false;
        break;
      }

      if (convertBondType(abp.getBondType()) != SINGLE) {
        //                if (charge_log && old_values[i] != 0)
        //                    fprintf(charge_log, "+%d*%c6", (abp->bond_type -
        //                    SINGLE), 'B' + (cip - ChargeIncTable));
        AtompKaValue[i] += (convertBondType(abp.getBondType()) - SINGLE) *
                           Options.ChargeIncTable[j].MultInc;
        //                if (charge_log && old_values[i] != 0)
        //                cip->mult_inc_used++;
      }

      // fetch alpha charge increment
      if (aap.getFormalCharge() != 0) {
        found = false;
        for (const auto &k : Options.ChargeIncTable)
          if (AtomSymbolMatch(aap.getSymbol(), k.AtomSymbol)) {
            //                if (charge_log && old_values[i] != 0)
            //                    fprintf(charge_log, "+%d*%c4", aap->charge,
            //                    'B' + (cip - ChargeIncTable));
            AtompKaValue[i] += aap.getFormalCharge() * k.AlphaInc;
            //                if (charge_log && old_values[i] != 0)
            //                cip->alpha_inc_used++;
            found = true;
            break;
          }
        if (!found) {
          // sprintf(msg_buffer, "%10s: no alpha increment for atom %d",
          // mp->name, i + 1);
          // AddMsgToList(msg_buffer);
          result = false;
          break;
        }
      }

      // fetch alpha acidity increment
      found = false;
      for (const auto &k : Options.AtomAcidity)
        if (AtomSymbolMatch(aap.getSymbol(), k.AtomSymbol)) {
          //            if (charge_log && old_values[i] != 0)
          //                fprintf(charge_log, "+%c16", 'B' + (ppa -
          //                alpha_path_table));
          //            if (charge_log && old_values[i] != 0)
          //                fprintf(charge_log, "*%c11", 'B' + (ip -
          //                atom_acidity_table));
          AtompKaValue[i] += ppa->Cond * k.AlphaInc;
          //            if (charge_log && old_values[i] != 0) ppa->cond_used++;
          //            if (charge_log && old_values[i] != 0)
          //            ip->alpha_inc_used++;
          found = true;
          break;
        }
      if (!found) {
        //                sprintf(msg_buffer, "%10s: no alpha increment for atom
        //                %d", mp->name, i + 1);
        //                AddMsgToList(msg_buffer);
        result = false;
        break;
      }

      // for all beta neighbours
      const Neighbourhood &nbph = neighbour_array[nbp.Atoms[j]];
      for (unsigned k = 0; k < nbph.Atoms.size(); k++) {
        if (nbph.Atoms[k] == i) continue;  // no loop back

        const Atom &bap = *Mol.getAtomWithIdx(nbph.Atoms[k]);  // beta atom
        const Bond &bbp = *Mol.getBondWithIdx(nbph.Bonds[k]);  // beta bond
        const PathEntry *ppb = nullptr;
        // fetch beta conductivity
        found = false;
        for (const auto &l : Options.BetaPathTable)  // ppb++
          if (AtomSymbolMatch(ap->getSymbol(), l.Path.AtomSymbol) &&
              LigandMatches(aap, abp, ppb->Path.Ligands[0], false) &&
              LigandMatches(bap, bbp, ppb->Path.Ligands[1], false)) {
            ppb = &Options.BetaPathTable[k];
            found = true;
            break;
          }
        if (!found) {
          //                    sprintf(msg_buffer, "%10s: no beta increment for
          //                    atom %d", mp->name, i + 1);
          //                    AddMsgToList(msg_buffer);
          result = false;
          break;
        }

        // fetch beta acidity increment
        found = false;
        for (const auto &l : Options.AtomAcidity)  // ip++
          if (AtomSymbolMatch(bap.getSymbol(), l.AtomSymbol)) {
            //                if (charge_log && old_values[i] != 0)
            //                    fprintf(charge_log, "+%c20", (int)('B' + (ppb
            //                    - beta_path_table)));
            //                if (charge_log && old_values[i] != 0)
            //                    fprintf(charge_log, "*%c12", 'B' + (ip -
            //                    atom_acidity_table));
            AtompKaValue[i] += ppb->Cond * l.BetaInc;
            //                if (charge_log && old_values[i] != 0)
            //                ppb->cond_used++;
            //                if (charge_log && old_values[i] != 0)
            //                ip->beta_inc_used++;
            found = true;
            break;
          }
        if (!found) {
          //                    sprintf(msg_buffer, "%10s: no beta increment for
          //                    atom %d", mp->name, i + 1);
          //                    AddMsgToList(msg_buffer);
          result = false;
          break;
        }

        // fetch beta charge increment
        if (bap.getFormalCharge() != 0) {
          found = false;
          for (const auto &l : Options.ChargeIncTable)
            if (AtomSymbolMatch(bap.getSymbol(), l.AtomSymbol)) {
              //                    if (charge_log && old_values[i] != 0)
              //                        fprintf(charge_log, "+%d*%c5",
              //                        bap->charge, 'B' + (cip -
              //                        ChargeIncTable));
              AtompKaValue[i] += bap.getFormalCharge() * l.BetaInc;  // ap
              found = true;
              break;
            }
          if (!found) {
            //                        sprintf(msg_buffer, "%10s: no beta
            //                        increment for atom %d", mp->name, i + 1);
            //                        AddMsgToList(msg_buffer);
            result = false;
            break;
          }
        }
      }
    }
    //        if (charge_log && old_values[i] != 0.0) fprintf(charge_log,
    //        "\t%g\n", ap->value);
    // The function is pKa = 7 + (pKa'-7)*beta + ((pKa' - 7)*alpha) ^ 3.
    AtompKaValue[i] = 7 + Options.Beta * (AtompKaValue[i] - 7) +
                      Options.Alpha * (AtompKaValue[i] - 7) * Options.Alpha *
                          (AtompKaValue[i] - 7) * Options.Alpha *
                          (AtompKaValue[i] - 7);
  }
  return result;
}

/*
 * Marks those atoms which have the smallest pKa value of all acidic
 * (i.e. already marked) atoms. All other marks are removed.
 * *pKa_value is set to the minimum pKa, and *gap is set to the
 * pKa difference to the next acidic set of atoms.
 */
int ChargeFix::markMostAcidicAtoms(double &pKa_value, double &gap) {
  double min_pKa = 1000.0;
  double next_pKa = 1000.0;
  int result = 0;
  double epsilon = 0.000001;
  unsigned na = Mol.getNumAtoms();

  for (unsigned i = 0; i < na; i++) {
    if (AtomColor[i] != 0 && AtompKaValue[i] < min_pKa)
      min_pKa = AtompKaValue[i];
  }
  for (unsigned i = 0; i < na; i++)
    if (AtomColor[i] != 0 && AtompKaValue[i] < min_pKa + epsilon) {
      result++;
      AtomColor[i] = 1;
    } else {
      if (AtomColor[i] != 0 && AtompKaValue[i] < next_pKa)
        next_pKa = AtompKaValue[i];
      AtomColor[i] = 0;
      AtompKaValue[i] = 0.;
    }

  pKa_value = min_pKa;
  gap = next_pKa - min_pKa;
  return result;
}

/*
 * Decrements the charges of all marked atoms.
 */
void ChargeFix::decrementMarkedCharges() {
  for (unsigned i = 0; i < Mol.getNumAtoms(); i++)
    if (0 != AtomColor[i])
      Mol.getAtomWithIdx(i)->setFormalCharge(
          Mol.getAtomWithIdx(i)->getFormalCharge() - 1);
}

struct atom_rank_t {
  unsigned atom;
  unsigned rank;
  unsigned n_ligands;
  double elneg;
  unsigned rank_sum;
};

/*
 * Refines the class of acidic atoms to a smaller one based on
 * the electronegativity of the neighbouring atoms.
 * It returns the size of this class.
 */
int ChargeFix::refineAcidicAtoms(std::vector<unsigned> &numbering) {
  int result = 0;
  bool changed, do_cis_trans;
  std::vector<atom_rank_t> atom_ranks(Mol.getNumAtoms());
  atom_rank_t ar_tmp;
  double epsilon = 0.000001;

  // set electronegativities
  for (unsigned i = 0; i < atom_ranks.size(); i++) {
    atom_ranks[i].atom = i + 1;
    atom_ranks[i].rank = 0;
    atom_ranks[i].rank_sum = 0;
    auto elneg = Options.ElnegTable.find(Mol.getAtomWithIdx(i)->getAtomicNum());
    if (Options.ElnegTable.end() != elneg)
      atom_ranks[i].elneg = elneg->second +
                            3.0 * Mol.getAtomWithIdx(i)->getFormalCharge() -
                            Options.Elneg0;  // elneg_table[0].value;
    else {
      // fprintf(stderr, "atom symbol '%s' not in periodic table\n",
      // mp->atom_array[i].AtomSymbol);
      return -1;
    }
  }
  // do preliminary ranking based on el. neg.
  for (unsigned i = 1; i < atom_ranks.size(); i++)  // sort by decreasing elneg.
    for (unsigned j = i; j > 0; j--)
      if (atom_ranks[j].elneg > atom_ranks[j - 1].elneg + epsilon) {
        ar_tmp = atom_ranks[j];
        atom_ranks[j] = atom_ranks[j - 1];
        atom_ranks[j - 1] = ar_tmp;
      } else
        break;
  atom_ranks[0].rank = 0;  // set ranks
  for (unsigned i = 1, j = 0; i < atom_ranks.size(); i++) {
    if (atom_ranks[i].elneg < atom_ranks[i - 1].elneg - epsilon) j = i;
    atom_ranks[i].rank = j;
  }
  for (unsigned i = 1; i < atom_ranks.size(); i++)  // resort by atom number
    for (unsigned j = i; j > 0; j--)
      if (atom_ranks[j].atom < atom_ranks[j - 1].atom) {
        ar_tmp = atom_ranks[j];
        atom_ranks[j] = atom_ranks[j - 1];
        atom_ranks[j - 1] = ar_tmp;
      } else
        break;

  // use unsaturation to split ranks (rank sum is misused here)
  for (auto &atom_rank : atom_ranks) atom_rank.rank_sum = 0;
  for (unsigned i = 0; i < Mol.getNumBonds(); i++) {
    const Bond &bond = *Mol.getBondWithIdx(i);
    if (convertBondType(bond.getBondType()) != SINGLE) {
      atom_ranks[bond.getBeginAtomIdx()].rank_sum++;
      atom_ranks[bond.getEndAtomIdx()].rank_sum++;
    }
  }
  for (unsigned i = 1; i < atom_ranks.size();
       i++)  // sort by rank(asc.) + unsat.(desc.)
    for (unsigned j = i; j > 0; j--)
      if (atom_ranks[j].rank < atom_ranks[j - 1].rank ||
          (atom_ranks[j].rank == atom_ranks[j - 1].rank &&
           atom_ranks[j].rank_sum > atom_ranks[j - 1].rank_sum)) {
        ar_tmp = atom_ranks[j];
        atom_ranks[j] = atom_ranks[j - 1];
        atom_ranks[j - 1] = ar_tmp;
      } else
        break;
  for (unsigned i = 1, j = 0; i < atom_ranks.size(); i++) {  // set new ranks
    if (atom_ranks[i].rank > atom_ranks[i - 1].rank ||
        atom_ranks[i].rank_sum < atom_ranks[i - 1].rank_sum)
      j = i;
    atom_ranks[i].rank = j;
  }
  for (unsigned i = 1; i < atom_ranks.size(); i++)  // restore atom number order
    for (unsigned j = i; j > 0; j--)
      if (atom_ranks[j].atom < atom_ranks[j - 1].atom) {
        ar_tmp = atom_ranks[j];
        atom_ranks[j] = atom_ranks[j - 1];
        atom_ranks[j - 1] = ar_tmp;
      } else
        break;

  for (auto &atom_rank : atom_ranks) atom_rank.n_ligands = 0;
  for (unsigned i = 0; i < Mol.getNumBonds(); i++) {
    const Bond &bond = *Mol.getBondWithIdx(i);
    atom_ranks[bond.getBeginAtomIdx()].n_ligands++;
    atom_ranks[bond.getEndAtomIdx()].n_ligands++;
  }

  // refine ranking using neighbour rank sums
  do_cis_trans = false;
  do {
    for (unsigned i = 0; i < atom_ranks.size(); i++)
      numbering[i] = atom_ranks[i].rank;
    // compute rank sums
    for (auto &atom_rank : atom_ranks) atom_rank.rank_sum = 0;
    for (unsigned i = 0; i < Mol.getNumBonds(); i++) {
      const Bond &bond = *Mol.getBondWithIdx(i);
      atom_ranks[bond.getBeginAtomIdx()].rank_sum +=
          atom_ranks[bond.getEndAtomIdx()].rank;
      atom_ranks[bond.getEndAtomIdx()].rank_sum +=
          atom_ranks[bond.getBeginAtomIdx()].rank;
    }
    if (do_cis_trans) {
      std::vector<RDGeom::Point3D> atomPoint;
      if (getMolAtomPoints(Mol, atomPoint))
        CisTransPerception(Mol, atomPoint, numbering, BondColor);

      for (unsigned i = 0; i < BondColor.size(); i++) {
        const Bond &bond = *Mol.getBondWithIdx(i);
        if (BondColor[i] != 0) {
          atom_ranks[bond.getBeginAtomIdx()].rank_sum += BondColor[i];
          atom_ranks[bond.getEndAtomIdx()].rank_sum += BondColor[i];
        }
      }
    }
    for (auto &atom_rank : atom_ranks)  // use average rank sum
      if (atom_rank.n_ligands > 0) {
        atom_rank.rank_sum *= 10;  // shift dec. point
        atom_rank.rank_sum /= atom_rank.n_ligands;
      }
    for (unsigned i = 1; i < atom_ranks.size(); i++)  // sort by rank + ranksum
      for (unsigned j = i; j > 0; j--)
        if (atom_ranks[j].rank < atom_ranks[j - 1].rank ||
            (atom_ranks[j].rank == atom_ranks[j - 1].rank &&
             atom_ranks[j].rank_sum < atom_ranks[j - 1].rank_sum)) {
          ar_tmp = atom_ranks[j];
          atom_ranks[j] = atom_ranks[j - 1];
          atom_ranks[j - 1] = ar_tmp;
        } else
          break;
    // compute new ranks
    changed = false;
    for (unsigned i = 1, j = 0; i < atom_ranks.size(); i++) {
      if (atom_ranks[i].rank > atom_ranks[i - 1].rank ||
          atom_ranks[i].rank_sum > atom_ranks[i - 1].rank_sum)
        j = i;
      changed = changed || atom_ranks[i].rank != j;
      atom_ranks[i].rank = j;
    }
    for (unsigned i = 1; i < atom_ranks.size();
         i++)  // restore atom number order
      for (unsigned j = i; j > 0; j--)
        if (atom_ranks[j].atom < atom_ranks[j - 1].atom) {
          ar_tmp = atom_ranks[j];
          atom_ranks[j] = atom_ranks[j - 1];
          atom_ranks[j - 1] = ar_tmp;
        } else
          break;
    if (!changed && !do_cis_trans) {
      do_cis_trans = true;
      changed = true;
    } else
      do_cis_trans = false;
  } while (changed);

  // find smalles rank of coloured atoms
  size_t min_rank = atom_ranks.size();
  for (unsigned i = 0; i < atom_ranks.size(); i++)
    if (AtomColor[i] != 0 && atom_ranks[i].rank < min_rank)
      min_rank = atom_ranks[i].rank;
  for (unsigned i = 0; i < atom_ranks.size(); i++) {
    if (AtomColor[i] != 0 &&
        atom_ranks[i].rank == min_rank) {  // count members of minimum class
      result++;
      AtompKaValue[i] = atom_ranks[i].rank + 1;
    } else {  // uncolour non-minimum members
      AtomColor[i] = 0;
      AtompKaValue[i] = 0.;
    }
  }
  if (result > 1) {
    for (unsigned int i : AtomColor)
      if (i != 0) {
        //                sprintf(msg_buffer, "atom %d in minimal rank class", i
        //                + 1);
        //                AddMsgToList(msg_buffer);
      }
  }
  return result;
}

}  // namespace StructureCheck
}  // namespace RDKit
