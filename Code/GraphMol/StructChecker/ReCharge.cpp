//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "StructChecker.h"
#include "Pattern.h"
#include "Utilites.h"
#include "ReCharge.h"

namespace RDKit {
 namespace StructureCheck {

/*
* Returns the total charge of all atoms in molecule.
*/
    int  TotalCharge     (const ROMol &mol) {
        int charge = 0;
        for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
            const Atom &atom = *mol.getAtomWithIdx(i);
            charge += atom.getFormalCharge();
        }
        return 0;
    }

/*
* Estimates the pKa value of acidic atoms in molecule *mp and sets
* the 'value' fields to this value.
* It returns TRUE if all went well and FALSE otherwise.
*/
bool SetpKaValues(const ROMol &mol, std::vector<double> &atom_pKa) {
    bool result = true;
    std::vector<Neighbourhood> neighbour_array(mol.getNumAtoms()); //, *nbp, *nbph;
    std::vector<unsigned> match;    //[MAXNEIGHBOURS + 1];

    AugmentedAtom *AAp;
    struct inc_entry_t *cip, *ip;
    struct path_entry_t *ppa, *ppb;

    struct reaccs_atom_t *ap, *aap, *bap;
    struct reaccs_bond_t *abp, *bbp;

// !!!
    static float old_values[MAXATOMS];


    SetupNeighbourhood(mp, neighbour_array, mp->n_atoms);

    ResetColors(mp);
    for (i = 0; i<mp->n_atoms; i++)
        old_values[i] = mp->atom_array[i].value;
    ResetValues(mp);

    /* for all atoms in *mp */
    for (i = 0, ap = mp->atom_array; i<mp->n_atoms; i++, ap++)
    {                                         /* Is atom acidic? */
        for (j = 0, AAp = acidic_atoms; j<nacidic; j++, AAp++)
            if (AAMatch(mp, i, match, AAp, (int *)NULL, neighbour_array))
                break;
        if (j == nacidic) continue; /* center not acidic -> goto next atom */
        ap->color++;
        ap->value = 0.0;

        if (charge_log && old_values[i] != 0.0)
            StartPredictionLine(AAp->short_name, i + 1, old_values[i]);
        if (old_values[i] != 0.0)
            fprintf(stderr, "atom %d: '%s'\n", i + 1, AAp->short_name);

        /* add local charge increment */
        for (j = 0, cip = charge_inc_table; j<ncharge; j++, cip++)
            if (AtomSymbolMatch(ap->atom_symbol, cip->atom_symbol))
            {
                ap->value += ap->charge*cip->local_inc;
                if (charge_log && old_values[i] != 0 && ap->charge != NONE) cip->local_inc_used++;
                if (charge_log && old_values[i] != 0) fprintf(charge_log, "+%d*%c3", ap->charge, 'B' + j);
                break;
            }

        /* add local atom acidity */
        for (j = 0, ip = atom_acidity_table; j<natomacidity; j++, ip++)
            if (AtomSymbolMatch(ap->atom_symbol, ip->atom_symbol))
            {
                ap->value += ip->local_inc;
                if (charge_log && old_values[i] != 0) ip->local_inc_used++;
                if (charge_log && old_values[i] != 0)
                    fprintf(charge_log, "+%c10", 'B' + j);
                break;
            }

        nbp = &neighbour_array[i];      /* for all alpha neighbours */
        for (j = 0; j<nbp->n_ligands; j++)
        {
            aap = &mp->atom_array[nbp->atoms[j]];  /* alpha atom */
            abp = &mp->bond_array[nbp->bonds[j]];  /* alpha bond */
                                                   /* fetch alpha path conductivity */
            for (k = 0, ppa = alpha_path_table; k<nalphapath; k++, ppa++)
                if (AtomSymbolMatch(ap->atom_symbol,
                    ppa->path.atom_symbol) &&
                    LigandMatches(aap, abp, &ppa->path.ligands[0], FALSE) &&
                    TRUE)
                    break;
            if (k == nalphapath)
            {
                sprintf(msg_buffer, "%10s: no alpha path for atom %d", mp->name, i + 1);
                AddMsgToList(msg_buffer);
                result = FALSE;
                break;
            }

            if (abp->bond_type != SINGLE)
            {
                if (charge_log && old_values[i] != 0)
                    fprintf(charge_log, "+%d*%c6",
                        (abp->bond_type - SINGLE),
                        'B' + (cip - charge_inc_table));
                ap->value += (abp->bond_type - SINGLE)*cip->mult_inc;
                if (charge_log && old_values[i] != 0) cip->mult_inc_used++;
            }

            /* fetch alpha charge increment */
            if (aap->charge != 0)
            {
                for (k = 0, cip = charge_inc_table;
                k<ncharge;
                    k++, cip++)
                    if (AtomSymbolMatch(aap->atom_symbol, cip->atom_symbol))
                        break;
                if (k == ncharge)
                {
                    sprintf(msg_buffer, "%10s: no alpha increment for atom %d", mp->name, i + 1);
                    AddMsgToList(msg_buffer);
                    result = FALSE;
                    break;
                }
                if (charge_log && old_values[i] != 0)
                    fprintf(charge_log, "+%d*%c4", aap->charge, 'B' + (cip - charge_inc_table));
                ap->value += aap->charge*cip->alpha_inc;
                if (charge_log && old_values[i] != 0) cip->alpha_inc_used++;
            }

            /* fetch alpha acidity increment */
            for (k = 0, ip = atom_acidity_table; k<natomacidity; k++, ip++)
                if (AtomSymbolMatch(aap->atom_symbol, ip->atom_symbol))
                    break;
            if (k == natomacidity)
            {
                sprintf(msg_buffer, "%10s: no alpha increment for atom %d", mp->name, i + 1);
                AddMsgToList(msg_buffer);
                result = FALSE;
                break;
            }
            if (charge_log && old_values[i] != 0)
                fprintf(charge_log, "+%c16", 'B' + (ppa - alpha_path_table));
            if (charge_log && old_values[i] != 0)
                fprintf(charge_log, "*%c11", 'B' + (ip - atom_acidity_table));
            ap->value += ppa->cond*ip->alpha_inc;
            if (charge_log && old_values[i] != 0) ppa->cond_used++;
            if (charge_log && old_values[i] != 0) ip->alpha_inc_used++;

            /* for all beta neighbours */
            nbph = &neighbour_array[nbp->atoms[j]];
            for (k = 0; k<nbph->n_ligands; k++)
            {
                if (nbph->atoms[k] == i) continue;  /* no loop back */

                bap = &mp->atom_array[nbph->atoms[k]]; /* beta atom */
                bbp = &mp->bond_array[nbph->bonds[k]]; /* beta bond */
                for (l = 0, ppb = beta_path_table; l<nbetapath; l++, ppb++) /* fetch beta conductivity */
                    if (AtomSymbolMatch(ap->atom_symbol, ppb->path.atom_symbol) &&
                        LigandMatches(aap, abp, &ppb->path.ligands[0], FALSE) &&
                        LigandMatches(bap, bbp, &ppb->path.ligands[1], FALSE) &&
                        TRUE)
                        break;
                if (l == nbetapath)
                {
                    sprintf(msg_buffer, "%10s: no beta increment for atom %d", mp->name, i + 1);
                    AddMsgToList(msg_buffer);
                    result = FALSE;
                    break;
                }

                /* fetch beta acidity increment */
                for (l = 0, ip = atom_acidity_table; l<natomacidity; l++, ip++)
                    if (AtomSymbolMatch(bap->atom_symbol, ip->atom_symbol))
                        break;
                if (l == natomacidity)
                {
                    sprintf(msg_buffer, "%10s: no beta increment for atom %d", mp->name, i + 1);
                    AddMsgToList(msg_buffer);
                    result = FALSE;
                    break;
                }
                if (charge_log && old_values[i] != 0)
                    fprintf(charge_log, "+%c20", (int)('B' + (ppb - beta_path_table)));
                if (charge_log && old_values[i] != 0)
                    fprintf(charge_log, "*%c12", 'B' + (ip - atom_acidity_table));
                ap->value += ppb->cond*ip->beta_inc;
                if (charge_log && old_values[i] != 0) ppb->cond_used++;
                if (charge_log && old_values[i] != 0) ip->beta_inc_used++;

                /* fetch beta charge increment */
                if (bap->charge != 0)
                {
                    for (l = 0, cip = charge_inc_table; l<ncharge; l++, cip++)
                        if (AtomSymbolMatch(bap->atom_symbol, cip->atom_symbol))
                            break;
                    if (l == ncharge)
                    {
                        sprintf(msg_buffer, "%10s: no beta increment for atom %d", mp->name, i + 1);
                        AddMsgToList(msg_buffer);
                        result = FALSE;
                        break;
                    }
                    if (charge_log && old_values[i] != 0)
                        fprintf(charge_log, "+%d*%c5", bap->charge, 'B' + (cip - charge_inc_table));
                    ap->value += bap->charge*cip->beta_inc;
                    if (charge_log && old_values[i] != 0) cip->beta_inc_used++;
                }
            }
        }
        if (charge_log && old_values[i] != 0.0) fprintf(charge_log, "\t%g\n", ap->value);
        ap->value = 7 + beta*(ap->value - 7)
            + alpha*(ap->value - 7)*
            alpha*(ap->value - 7)*
            alpha*(ap->value - 7);
    }

    MyFree((char *)neighbour_array);
    return (result);
}

int MarkMostAcidicAtoms(struct reaccs_molecule_t *mp,
    double *pKa_value, double *gap)
    /*
    * Marks those atoms which have the smallest pKa value of all acidic
    * (i.e. already marked) atoms. All other marks are removed.
    * *pKa_value is set to the minimum pKa, and *gap is set to the
    * pKa difference to the next acidic set of atoms.
    */
{
    double min_pKa, next_pKa;
    int i, result = 0;
    struct reaccs_atom_t *ap;
    double epsilon = 0.000001;

    for (i = 0, min_pKa = 1000.0, ap = mp->atom_array; i<mp->n_atoms; i++, ap++)
        if (ap->color != NONE  &&  ap->value < min_pKa)
            min_pKa = ap->value;

    next_pKa = 1000.0;
    for (i = 0, ap = mp->atom_array; i<mp->n_atoms; i++, ap++)
        if (ap->color != NONE  &&  ap->value < min_pKa + epsilon)
        {
            result++;
            ap->color = NONE + 1;
        }
        else
        {
            if (ap->color != NONE  &&  ap->value < next_pKa)
                next_pKa = ap->value;
            ap->color = NONE;
            ap->value = 0;
        }

    *pKa_value = min_pKa;
    *gap = next_pKa - min_pKa;
    return (result);
}

void DecrementMarkedCharges(struct reaccs_molecule_t *mp)
/*
* Decrements the charges of all marked atoms in *mp.
*/
{
    int i;

    for (i = 0; i<mp->n_atoms; i++)
        if (mp->atom_array[i].color != NONE)
            mp->atom_array[i].charge--;
}

struct atom_rank_t
{
    int atom;
    int rank;
    int n_ligands;
    double elneg;
    int rank_sum;
};

int RefineAcidicAtoms(struct reaccs_molecule_t *mp,
    int numbering[])
    /*
    * Refines the class of acidic atoms to a smaller one based on
    * the electronegativity of the neighbouring atoms.
    * It returns the size of this class.
    */
{
    int result = 0;
    int i, j, min_rank;
    int changed, do_cis_trans;
    struct atom_rank_t *atom_ranks, ar_tmp;
    struct reaccs_bond_t *bp;
    double epsilon = 0.000001;

    /* set electronegativities */
    atom_ranks = TypeAlloc(mp->n_atoms, struct atom_rank_t);
    for (i = 0; i<mp->n_atoms; i++)
    {
        atom_ranks[i].atom = i + 1;
        atom_ranks[i].rank = 0;
        atom_ranks[i].rank_sum = 0;
        for (j = 0; j<nelneg; j++)
            if (0 == strcmp(mp->atom_array[i].atom_symbol, elneg_table[j].symbol))
            {
                atom_ranks[i].elneg = elneg_table[j].value;
                break;
            }
        if (j == nelneg)
        {
            fprintf(stderr, "atom symbol '%s' not in periodic table\n", mp->atom_array[i].atom_symbol);
            MyFree((char *)atom_ranks);
            return (-1);
        }
        atom_ranks[i].elneg += 3.0*mp->atom_array[i].charge -
            elneg_table[0].value;
    }

    /* do preliminary ranking based on el. neg. */
    for (i = 1; i<mp->n_atoms; i++)    /* sort by decreasing elneg. */
        for (j = i; j>0; j--)
            if (atom_ranks[j].elneg > atom_ranks[j - 1].elneg + epsilon)
            {
                ar_tmp = atom_ranks[j];
                atom_ranks[j] = atom_ranks[j - 1];
                atom_ranks[j - 1] = ar_tmp;
            }
            else
                break;
    atom_ranks[0].rank = 0;            /* set ranks */
    for (i = 1, j = 0; i<mp->n_atoms; i++)
    {
        if (atom_ranks[i].elneg < atom_ranks[i - 1].elneg - epsilon) j = i;
        atom_ranks[i].rank = j;
    }
    for (i = 1; i<mp->n_atoms; i++)  /* resort by atom number */
        for (j = i; j>0; j--)
            if (atom_ranks[j].atom < atom_ranks[j - 1].atom)
            {
                ar_tmp = atom_ranks[j];
                atom_ranks[j] = atom_ranks[j - 1];
                atom_ranks[j - 1] = ar_tmp;
            }
            else
                break;

    /* use unsaturation to split ranks (rank sum is misused here) */
    for (i = 0; i<mp->n_atoms; i++) atom_ranks[i].rank_sum = 0;
    for (i = 0, bp = mp->bond_array; i<mp->n_bonds; i++, bp++)
        if (bp->bond_type != SINGLE)
        {
            atom_ranks[bp->atoms[0] - 1].rank_sum++;
            atom_ranks[bp->atoms[1] - 1].rank_sum++;
        }
    for (i = 1; i<mp->n_atoms; i++) /* sort by rank(asc.) + unsat.(desc.) */
        for (j = i; j>0; j--)
            if (atom_ranks[j].rank < atom_ranks[j - 1].rank ||
                (atom_ranks[j].rank == atom_ranks[j - 1].rank  &&
                    atom_ranks[j].rank_sum > atom_ranks[j - 1].rank_sum))
            {
                ar_tmp = atom_ranks[j];
                atom_ranks[j] = atom_ranks[j - 1];
                atom_ranks[j - 1] = ar_tmp;
            }
            else
                break;
    for (i = 1, j = 0; i<mp->n_atoms; i++) /* set new ranks */
    {
        if (atom_ranks[i].rank > atom_ranks[i - 1].rank ||
            atom_ranks[i].rank_sum < atom_ranks[i - 1].rank_sum) j = i;
        atom_ranks[i].rank = j;
    }
    for (i = 1; i<mp->n_atoms; i++)    /* restore atom number order */
        for (j = i; j>0; j--)
            if (atom_ranks[j].atom < atom_ranks[j - 1].atom)
            {
                ar_tmp = atom_ranks[j];
                atom_ranks[j] = atom_ranks[j - 1];
                atom_ranks[j - 1] = ar_tmp;
            }
            else
                break;

    for (i = 0; i<mp->n_atoms; i++) atom_ranks[i].n_ligands = 0;
    for (i = 0, bp = mp->bond_array; i<mp->n_bonds; i++, bp++)
    {
        atom_ranks[bp->atoms[0] - 1].n_ligands++;
        atom_ranks[bp->atoms[1] - 1].n_ligands++;
    }

    /* refine ranking using neighbour rank sums */
    do_cis_trans = FALSE;
    do
    {
        for (i = 0; i<mp->n_atoms; i++)
            numbering[i] = atom_ranks[i].rank;
        /* compute rank sums */
        for (i = 0; i<mp->n_atoms; i++) atom_ranks[i].rank_sum = 0;
        for (i = 0, bp = mp->bond_array; i<mp->n_bonds; i++, bp++)
        {
            atom_ranks[bp->atoms[0] - 1].rank_sum +=
                atom_ranks[bp->atoms[1] - 1].rank;
            atom_ranks[bp->atoms[1] - 1].rank_sum +=
                atom_ranks[bp->atoms[0] - 1].rank;
        }
        if (do_cis_trans)
        {
            CisTransPerception(mp, numbering);
            for (i = 0, bp = mp->bond_array; i<mp->n_bonds; i++, bp++)
                if (bp->color != NONE)
                {
                    atom_ranks[bp->atoms[0] - 1].rank_sum += bp->color;
                    atom_ranks[bp->atoms[1] - 1].rank_sum += bp->color;
                }
        }
        for (i = 0; i<mp->n_atoms; i++) /* use average rank sum */
            if (atom_ranks[i].n_ligands > 0)
            {
                atom_ranks[i].rank_sum *= 10; /* shift dec. point */
                atom_ranks[i].rank_sum /= atom_ranks[i].n_ligands;
            }
        for (i = 1; i<mp->n_atoms; i++)   /* sort by rank + ranksum */
            for (j = i; j>0; j--)
                if (atom_ranks[j].rank < atom_ranks[j - 1].rank ||
                    (atom_ranks[j].rank == atom_ranks[j - 1].rank  &&
                        atom_ranks[j].rank_sum < atom_ranks[j - 1].rank_sum))
                {
                    ar_tmp = atom_ranks[j];
                    atom_ranks[j] = atom_ranks[j - 1];
                    atom_ranks[j - 1] = ar_tmp;
                }
                else
                    break;
        /* compute new ranks */
        for (i = 1, changed = FALSE, j = 0; i<mp->n_atoms; i++)
        {
            if (atom_ranks[i].rank > atom_ranks[i - 1].rank ||
                atom_ranks[i].rank_sum > atom_ranks[i - 1].rank_sum) j = i;
            changed = changed || atom_ranks[i].rank != j;
            atom_ranks[i].rank = j;
        }
        for (i = 1; i<mp->n_atoms; i++) /* restore atom number order */
            for (j = i; j>0; j--)
                if (atom_ranks[j].atom < atom_ranks[j - 1].atom)
                {
                    ar_tmp = atom_ranks[j];
                    atom_ranks[j] = atom_ranks[j - 1];
                    atom_ranks[j - 1] = ar_tmp;
                }
                else
                    break;
        if (!changed && !do_cis_trans)
        {
            do_cis_trans = TRUE;
            changed = TRUE;
        }
        else
            do_cis_trans = FALSE;
    } while (changed);

    /* find smalles rank of coloured atoms */
    for (i = 0, min_rank = mp->n_atoms; i<mp->n_atoms; i++)
        if (mp->atom_array[i].color != NONE  &&
            atom_ranks[i].rank < min_rank)
            min_rank = atom_ranks[i].rank;
    for (i = 0; i<mp->n_atoms; i++)
    {
        if (mp->atom_array[i].color != NONE  &&
            atom_ranks[i].rank == min_rank)
        {                        /* count members of minimum class */
            result++;
            mp->atom_array[i].value = atom_ranks[i].rank + 1;
        }
        else
        {                        /* uncolour non-minimum members */
            mp->atom_array[i].color = NONE;
            mp->atom_array[i].value = 0;
        }
    }

    MyFree((char *)atom_ranks);
    if (result > 1)
    {
        for (i = 0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color != NONE)
            {
                sprintf(msg_buffer, "atom %d in minimal rank class", i + 1);
                AddMsgToList(msg_buffer);
            }
    }
    return (result);
}

/*
* Removes hydrogens from *mp until desired_charge is reached. The
* positions for hydrogen removal are selected by "acidity" combined
* with a refinement algorithm. It returns TRUE if molecule could be
* neutralized and FALSE if any problem were encountered.
* *ndeprot and *nrefine are set to the number of deprotonations
* and refinement cycles performed.
*/
    bool RechargeMolecule(RWMol &mol, int desired_charge, unsigned &ndeprot, unsigned &nrefine) {
        bool changed = false;
        int nacid;
        double gap, pKa_value;
        std::vector<double> atom_pKa (mol.getNumAtoms());
        double acidity_limit = 24.0;
        std::vector<unsigned> numbering(mol.getNumAtoms());

        ndeprot = 0;    // number of deprotonation cycles
        nrefine = 0;    // number of refinements necessary
        for (;;)
        {
            while (TotalCharge(mol) > desired_charge)
            {
                SetpKaValues(mol, atom_pKa);
                nacid = MarkMostAcidicAtoms(mol, pKa_value, gap, acidity_limit);
                if (pKa_value > acidity_limit) {
//                    sprintf(msg_buffer, "pKa_value = %.2g", pKa_value);
//                    AddMsgToList(msg_buffer);
                    return false;     // not acidic enough
                }
                else if (nacid == 1 || (nacid == TotalCharge(mol) && gap > 8.0)) { // easy case
//                    sprintf(msg_buffer, "pKa = %.2g", pKa_value);
//                    AddMsgToList(msg_buffer);
                    DecrementMarkedCharges(mol);
                }
                else
                    break;
                ndeprot++;
            }

            if (TotalCharge(mol) > desired_charge) {
                nrefine++;
//// ?? expl. H ??                numbering = TypeAlloc(mp->n_atoms, int);
                nacid = RefineAcidicAtoms(mol, numbering);
                if (nacid == 1)
                    DecrementMarkedCharges(mol);
                else if (AllCentersRefined(mol, numbering))
                {
                    for (unsigned i = 0; i < mol.getNumAtoms(); i++) // Select one mark
                        if (mp->atom_array[i].color != NONE && nacid != 1)
                        {
                            nacid--; mp->atom_array[i].color = NONE;
                        }
                    DecrementMarkedCharges(mol);
                }
                else
                {
//                    sprintf(msg_buffer, "%10s: could not fix charges", mp->name);
//                    AddMsgToList(msg_buffer);
                    return false;
                }
                ndeprot++;
            }
            else
            {
                ResetValues(mol);
                return true;
            }
        }
        return changed;
    }

 }// namespace StructureCheck
} // namespace RDKit
