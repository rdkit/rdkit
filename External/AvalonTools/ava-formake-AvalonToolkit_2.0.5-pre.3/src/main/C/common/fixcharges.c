//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

/************************************************************************/
/*                                                                      */
/*    File:           fixcharges.c                                      */
/*                                                                      */
/*    Purpose:        Module for fixing charges in molecules            */
/*                    which have a net charge. It updates the charges   */
/*                    of the acidic atoms such that a desired charge    */
/*                    given as a parameter is achieved.                 */
/*                                                                      */
/*                    The actual main module is for test purposes       */
/*                    only. The actual module exports only a few        */
/*                    functions doing the actual work.                  */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "forio.h"
#include "local.h"
#include "pattern.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "stereo.h"
#include "utilities.h"

#include "fixcharges.h"

static int SkipCommentLines(Fortran_FILE *fp)
{
   while (fp->status == FORTRAN_NORMAL  && fp->buffer[0] == '#')
      GetBuffer(fp);
   return (fp->status);
}

/*
 * Data structure to store the augmented atom descriptions for those
 * atoms considered acidic.
 */
#define MAXACIDIC     45
// static */
augmented_atom_t acidic_atoms[MAXACIDIC];
static int nacidic = 0;

static int ReadAcidicAtoms(Fortran_FILE *fp)
/*
 * Fills the acidic atom descriptions from file *fp into the local
 * data structure acidic_atoms[]. It returns TRUE if successfull and
 * FALSE if otherwise.
 */
{
   int i, j;
   char *aastart;
   int aalength;
   char buffer[255];
   augmented_atom_t *aap;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

   nacidic  = 0;
   sscanf(fp->buffer,"%d",&nacidic); GetBuffer(fp);

   if (nacidic < 1) return (FALSE);   /* no acidic atoms -> nonsense */
   if (nacidic > MAXACIDIC)
   {
      fprintf(stderr,"too many acidic atom\nincrease MAXMAXACIDIC\n");
      exit (EXIT_FAILURE);
   }

   for (i=0; i<nacidic  &&  fp->status == FORTRAN_NORMAL; i++)
   {
      if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (0);
      aastart = strchr(fp->buffer,'"')+1;
      if (!aastart || !strchr(aastart,'"')) return (0);
      aalength = (int)(strchr(aastart,'"') - aastart);
      strncpy(buffer,aastart,aalength); buffer[aalength] = '\0';
      if (!StringToAugmentedAtom(buffer,&acidic_atoms[i]))
         return(FALSE);
      aap = &acidic_atoms[i];
      for (j=0; j<aap->n_ligands; j++)
         if (aap->ligands[j].charge == NONE)
            aap->ligands[j].charge = ANY_CHARGE;
      GetBuffer(fp);
   }
   return (TRUE);
}

/*
 * Data structures used to store the local, alpha, and beta charge
 * increments. Those increments can be made dependent on the atom type.
 */
#define MAXCHARGEINC    5
/* static */ struct inc_entry_t
   {
      char atom_symbol[20];
      double local_inc;     int local_inc_used;
      double alpha_inc;     int alpha_inc_used;
      double beta_inc;      int beta_inc_used;
      double mult_inc;      int mult_inc_used;
   } charge_inc_table[MAXCHARGEINC];
static int ncharge = 0;

static int ReadChargeIncrements(Fortran_FILE *fp)
/*
 * Reads the local, alpha, and beta charge increments from file *fp
 * into the local data structure charge_inc_table[].
 * It returns TRUE on success and FALSE otherwise.
 */
{
   int i;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);
   if (sscanf(fp->buffer,"%d",&ncharge) != 1 ||
       ncharge < 0)
      return (FALSE);
   else
      GetBuffer(fp);

   for (i=0; i<ncharge  &&  SkipCommentLines(fp) == FORTRAN_NORMAL; i++)
   {
      if (sscanf(fp->buffer,"%s %lf %lf %lf %lf",
                            charge_inc_table[i].atom_symbol,
                            &charge_inc_table[i].local_inc,
                            &charge_inc_table[i].alpha_inc,
                            &charge_inc_table[i].beta_inc,
                            &charge_inc_table[i].mult_inc) != 5)
      return (FALSE);
      charge_inc_table[i].local_inc_used = 0;
      charge_inc_table[i].alpha_inc_used = 0;
      charge_inc_table[i].beta_inc_used = 0;
      charge_inc_table[i].mult_inc_used = 0;
      GetBuffer(fp);
   }
   if (fp->status == FORTRAN_NORMAL  ||  fp->status == FORTRAN_EOF)
      return (TRUE);
   else
      return (FALSE);
}

/*
 * Data structures to store atom dependent acidity factors.
 */
#define MAXATOMACIDITY 15
struct inc_entry_t atom_acidity_table[MAXATOMACIDITY];
int natomacidity = 0;

static int ReadAtomAcidityIncrements(Fortran_FILE *fp)
/*
 * Reads the local, alpha, and beta atom acidity increments from file *fp
 * into the local data structure atom_acidity_table[].
 * It returns TRUE on success and FALSE otherwise.
 */
{
   int i;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);
   if (sscanf(fp->buffer,"%d",&natomacidity) != 1 ||
       natomacidity < 0)
      return (FALSE);
   else
      GetBuffer(fp);

   for (i=0; i<natomacidity  &&  SkipCommentLines(fp) == FORTRAN_NORMAL; i++)
   {
      if (sscanf(fp->buffer,"%s %lf %lf %lf",
                            atom_acidity_table[i].atom_symbol,
                           &atom_acidity_table[i].local_inc,
                           &atom_acidity_table[i].alpha_inc,
                           &atom_acidity_table[i].beta_inc) != 4)
         return (FALSE);
      atom_acidity_table[i].local_inc_used = 0;
      atom_acidity_table[i].alpha_inc_used = 0;
      atom_acidity_table[i].beta_inc_used = 0;
      atom_acidity_table[i].mult_inc_used = 0;
      GetBuffer(fp);
   }
   if (fp->status == FORTRAN_NORMAL  ||  fp->status == FORTRAN_EOF)
      return (TRUE);
   else
      return (FALSE);
}

/*
 * Data structures to store alpha path conductivities. The paths are
 * stored in augmented atoms, which is a slight misuse.
 */
#define MAXALPHAPATH 10
// static */
struct path_entry_t
   {
      augmented_atom_t path;
      double cond;  int cond_used;
   } alpha_path_table[MAXALPHAPATH];
static int nalphapath = 0;

static int ReadAlphaPathTable(Fortran_FILE *fp)
/*
 * Reads the alpha charge conductivities from file *fp into the
 * local data structure alpha_path_table[]. The paths are stored
 * in an augmented atom which is a slight misuse of the structure.
 * The function returns TRUE on success and FALSE otherwise.
 */
{
   int i, j;
   char *aastart;
   int aalength;
   char buffer[255];
   augmented_atom_t *aap;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

   nalphapath = -1;
   sscanf(fp->buffer,"%d",&nalphapath); GetBuffer(fp);

   if (nalphapath < 0) return (FALSE);
   if (nalphapath > MAXALPHAPATH)
   {
      fprintf(stderr,"too many alpha path conductivities\n");
      fprintf(stderr,"increase MAXALPHAPATH\n");
      exit (EXIT_FAILURE);
   }
   for (i=0; i<nalphapath  &&  fp->status == FORTRAN_NORMAL; i++)
   {
      if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

      aastart = strchr(fp->buffer,'"')+1;   /* prepare augmented atom */
      if (!aastart || !strchr(aastart,'"')) return (0);
      aalength = (int)(strchr(aastart,'"') - aastart);
      strncpy(buffer,aastart,aalength); buffer[aalength] = '\0';

      if (1 != sscanf(aastart+aalength+1,"%lf",&alpha_path_table[i].cond)  ||
          !StringToAugmentedAtom(buffer,&alpha_path_table[i].path))
         return(FALSE);

      alpha_path_table[i].cond_used = 0;

      aap = &alpha_path_table[i].path;
      for (j=0; j<aap->n_ligands; j++)
         if (aap->ligands[j].charge == NONE)
            aap->ligands[j].charge = ANY_CHARGE;
      GetBuffer(fp);
   }
   return (TRUE);
}

/*
 * Data structures to store beta path conductivities. The paths are
 * stored in augmented atoms, which is a slight misuse.
 */
#define MAXBETAPATH    20
// static */
struct path_entry_t beta_path_table[MAXBETAPATH];
static int nbetapath = 0;

static int ReadBetaPathTable(Fortran_FILE *fp)
/*
 * Reads the beta charge conductivities from file *fp into the
 * local data structure beta_path_table[]. The paths are stored
 * in an augmented atom which is a slight misuse of the structure.
 * The function returns TRUE on success and FALSE otherwise.
 */
{
   int i, j;
   char *aastart;
   int aalength;
   char buffer[255];
   augmented_atom_t *aap;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

   nbetapath = -1;
   sscanf(fp->buffer,"%d",&nbetapath); GetBuffer(fp);

   if (nbetapath < 0) return (FALSE);
   if (nbetapath > MAXBETAPATH)
   {
      fprintf(stderr,"too many beta path conductivities\n");
      fprintf(stderr,"increase MAXBETAPATH\n");
      exit (EXIT_FAILURE);
   }
   for (i=0; i<nbetapath  &&  fp->status == FORTRAN_NORMAL; i++)
   {
      if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

      aastart = strchr(fp->buffer,'"')+1;   /* prepare augmented atom */
      if (!aastart || !strchr(aastart,'"')) return (0);
      aalength = (int)(strchr(aastart,'"') - aastart);
      strncpy(buffer,aastart,aalength); buffer[aalength] = '\0';

      if (1 != sscanf(aastart+aalength+1,"%lf",&beta_path_table[i].cond)  ||
          !StringToAugmentedAtom(buffer,&beta_path_table[i].path))
         return(FALSE);

      beta_path_table[i].cond_used = 0;

      aap = &beta_path_table[i].path;
      for (j=0; j<aap->n_ligands; j++)
      if (aap->ligands[j].charge == NONE)
         aap->ligands[j].charge = ANY_CHARGE;
      GetBuffer(fp);
   }
   return (TRUE);
}

static double alpha, beta;

static int ReadTransformation(Fortran_FILE *fp)
/*
 * Reads the coefficients of the transfomation function used
 * to stretch the preliminary pKa values to the actual predictions.
 * The function is pKa = 7 + (pKa'-7)*beta + ((pKa'-7)*alpha)^3.
 */
{
   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (FALSE);

   if (2 != sscanf(fp->buffer,"%lf %lf",&alpha, &beta))
      return (FALSE);
   else
   {
      GetBuffer(fp);
      return (TRUE);
   }
}

// static */
struct elneg_entry_t
   {
      char symbol[4];
      double value;
   } elneg_table[110];
static int nelneg = 0;

static int ReadElNegTable(Fortran_FILE *fp, struct elneg_entry_t elneg_table[])
/*
 * Reads a table of electronegativities consisting of rows with pairs
 * chemical symbol vs. el. neg.
 */
{
   int nelneg;
   int i;

   if (SkipCommentLines(fp) != FORTRAN_NORMAL) return (0);

   if (sscanf(fp->buffer,"%d",&nelneg) != 1) return (0);
   else                                      GetBuffer(fp);

   for (i=0; i<nelneg  &&  SkipCommentLines(fp) == FORTRAN_NORMAL; i++)
   {
      if (sscanf(fp->buffer,"%s %lf",elneg_table[i].symbol,
                                    &elneg_table[i].value) != 2)
         return (i);
      GetBuffer(fp);
   }
   if (fp->status == FORTRAN_NORMAL  ||
       fp->status == FORTRAN_EOF)
      return (nelneg);
   else
      return (0);
}

char pKa_version[10] = {'\0'};

int InitializeChargeDataTables(Fortran_FILE *fp)
/*
 * Initializes the own variables of the charge computation
 * module from file *fp. Returns TRUE if everything went ok
 * and FALSE otherwise.
 */
{
   char *cp, *cph;

   cp = strchr(fp->buffer, '_');
   if (cp)
   {
      cph = strchr(cp+1, '_');
      if (cph)
      {
         strncpy(pKa_version, cp+1, (int)(cph-cp)-1);
         pKa_version[(int)(cph-cp)-1] = '\0';
      }
      else
      {
         strncpy(pKa_version, cp+1, 9);
         pKa_version[9] = '\0';
      }
      if (log_file) fprintf(stderr,"pKa version = %s\n",pKa_version);
   }

   if (!ReadAcidicAtoms(fp))
   {
      fprintf(stderr, "Error reading acidic atoms\n");
      exit (EXIT_FAILURE);
   }

   if (!ReadChargeIncrements(fp))
   {
      fprintf(stderr, "Error reading charge increments\n");
      exit (EXIT_FAILURE);
   }

   if (!ReadAtomAcidityIncrements(fp))
   {
      fprintf(stderr, "Error reading atom acidity increments\n");
      exit (EXIT_FAILURE);
   }

   if (!ReadAlphaPathTable(fp))
   {
      fprintf(stderr, "Error reading alpha path table\n");
      exit (EXIT_FAILURE);
   }

   if (!ReadBetaPathTable(fp))
   {
      fprintf(stderr, "Error reading beta path table\n");
      exit (EXIT_FAILURE);
   }

   if (!ReadTransformation(fp))
   {
      fprintf(stderr, "Error reading transformation polynom\n");
      exit (EXIT_FAILURE);
   }
   nelneg = ReadElNegTable(fp, elneg_table);

   if (nelneg == 0)
      return (FALSE);
   else
      return (TRUE);
}

/*
 * The following routines are used to support the fitting of charge
 * parameters. The idea is to write a tab separated table of
 * multipliers for optimisation using the Excel solver facility.
 * Each multiplier line is prepended by a message and the atom number.
 */

static FILE *charge_log=NULL;
void SetChargeLog(FILE *logfile)
{
   charge_log = logfile;
}

static char charge_message[255];
void SetChargeMessage(char *message)
{
   strncpy(charge_message, message, 254);
}

void StartPredictionLine(char *aa_string, int atom, double value)
/*
 * Prints the header of an atom prediction line.
 */
{
   if (!charge_log) return;
   fprintf(charge_log, "%s\t%d\t%s\t%g\t=0",
                       charge_message,atom,aa_string,value);
}

void PrintChargeHeader()
{
   int i;

   if (!charge_log) return;

   fprintf(charge_log,"Charge Increments\n");
   fprintf(charge_log,"symbol");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%s",charge_inc_table[i].atom_symbol);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"local_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%g",charge_inc_table[i].local_inc);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"alpha_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%g",charge_inc_table[i].alpha_inc);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"beta_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%g",charge_inc_table[i].beta_inc);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"mult_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%g",charge_inc_table[i].mult_inc);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Atom Acidity Factors\n");
   fprintf(charge_log,"symbol");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%s",atom_acidity_table[i].atom_symbol);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"local_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%g",atom_acidity_table[i].local_inc);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"alpha_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%g",atom_acidity_table[i].alpha_inc);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"beta_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%g",atom_acidity_table[i].beta_inc);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Alpha Path Conductivity\n");
   fprintf(charge_log,"Path");
   for (i=0; i<nalphapath; i++)
      fprintf(charge_log,"\t%s",alpha_path_table[i].path.short_name);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"conductivity");
   for (i=0; i<nalphapath; i++)
      fprintf(charge_log,"\t%g",alpha_path_table[i].cond);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Beta Path Conductivity\n");
   fprintf(charge_log,"Path");
   for (i=0; i<nbetapath; i++)
      fprintf(charge_log,"\t%s",beta_path_table[i].path.short_name);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"conductivity");
   for (i=0; i<nbetapath; i++)
      fprintf(charge_log,"\t%g",beta_path_table[i].cond);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log, "Molecule\tCenter\tAugmented Atom\tpKa\tComputation\tPrediction\n");
}

void PrintChargeFooter()
{
   int i;

   if (!charge_log) return;

   fprintf(charge_log,"\n\nCharge Increments\n");
   fprintf(charge_log,"symbol");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%s",charge_inc_table[i].atom_symbol);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"local_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%d",charge_inc_table[i].local_inc_used);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"alpha_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%d",charge_inc_table[i].alpha_inc_used);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"beta_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%d",charge_inc_table[i].beta_inc_used);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"mult_inc");
   for (i=0; i<ncharge; i++)
      fprintf(charge_log,"\t%d",charge_inc_table[i].mult_inc_used);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Atom Acidity Factors\n");
   fprintf(charge_log,"symbol");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%s",atom_acidity_table[i].atom_symbol);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"local_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%d",atom_acidity_table[i].local_inc_used);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"alpha_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%d",atom_acidity_table[i].alpha_inc_used);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"beta_inc");
   for (i=0; i<natomacidity; i++)
      fprintf(charge_log,"\t%d",atom_acidity_table[i].beta_inc_used);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Alpha Path Conductivity\n");
   fprintf(charge_log,"Path");
   for (i=0; i<nalphapath; i++)
      fprintf(charge_log,"\t%s",alpha_path_table[i].path.short_name);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"conductivity");
   for (i=0; i<nalphapath; i++)
      fprintf(charge_log,"\t%d",alpha_path_table[i].cond_used);
   fprintf(charge_log,"\n\n");

   fprintf(charge_log,"Beta Path Conductivity\n");
   fprintf(charge_log,"Path");
   for (i=0; i<nbetapath; i++)
      fprintf(charge_log,"\t%s",beta_path_table[i].path.short_name);
   fprintf(charge_log,"\n");
   fprintf(charge_log,"conductivity");
   for (i=0; i<nbetapath; i++)
      fprintf(charge_log,"\t%d",beta_path_table[i].cond_used);
   fprintf(charge_log,"\n\n");
}

static int LigandMatches(struct reaccs_atom_t *ap,
                         struct reaccs_bond_t *bp,
                         ligand_t             *lp,
                         int                   use_charge)
/*
 * Checks if the ligand *lp is compatible with the bond *bp
 * and the atom *ap.
 */
{
   if (lp->bond_type != ANY_BOND  &&
       bp->bond_type != lp->bond_type) return (FALSE);

   if (lp->radical != ANY_RADICAL  &&
       ap->radical != lp->radical) return (FALSE);
   if ((lp->charge != NONE || use_charge) &&
       lp->charge != ANY_CHARGE       &&
       ap->charge != lp->charge) return (FALSE);

   return (AtomSymbolMatch(ap->atom_symbol, lp->atom_symbol));
}

int SetpKaValues(struct reaccs_molecule_t *mp)
/*
 * Estimates the pKa value of acidic atoms in molecule *mp and sets
 * the 'value' fields to this value.
 * It returns TRUE if all went well and FALSE otherwise.
 */
{
   int result;
   int i, j, k, l;
   neighbourhood_t *neighbour_array, *nbp, *nbph;
   unsigned int match[MAXNEIGHBOURS+1];

   augmented_atom_t *AAp;
   struct inc_entry_t *cip, *ip;
   struct path_entry_t *ppa, *ppb;

   struct reaccs_atom_t *ap, *aap, *bap;
   struct reaccs_bond_t *abp, *bbp;

   static float old_values[MAXATOMS];

   result = TRUE;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   ResetColors(mp);
   for (i=0; i<mp->n_atoms; i++)
     old_values[i] = mp->atom_array[i].value;
   ResetValues(mp);

                                            /* for all atoms in *mp */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {                                         /* Is atom acidic? */
      for (j=0, AAp=acidic_atoms; j<nacidic; j++, AAp++)
         if (AAMatch(mp,i,match,AAp,(int *)NULL,neighbour_array))
            break;
      if (j == nacidic) continue; /* center not acidic -> goto next atom */
      ap->color++;
      ap->value = 0.0;

      if (charge_log && old_values[i] != 0.0)
         StartPredictionLine(AAp->short_name, i+1, old_values[i]);
//      if (old_values[i] != 0.0)
//         fprintf(stderr,"atom %d: '%s'\n", i+1, AAp->short_name);

      /* add local charge increment */
      for (j=0, cip=charge_inc_table; j<ncharge; j++,cip++)
         if (AtomSymbolMatch(ap->atom_symbol, cip->atom_symbol))
         {
            ap->value += ap->charge*cip->local_inc;
            if (charge_log && old_values[i] != 0 && ap->charge != NONE) cip->local_inc_used++;
            if (charge_log && old_values[i] != 0) fprintf(charge_log,"+%d*%c3",ap->charge,'B'+j);
            break;
         }

                               /* add local atom acidity */
      for (j=0, ip=atom_acidity_table; j<natomacidity; j++,ip++)
         if (AtomSymbolMatch(ap->atom_symbol, ip->atom_symbol))
         {
            ap->value += ip->local_inc;
            if (charge_log && old_values[i] != 0) ip->local_inc_used++;
            if (charge_log && old_values[i] != 0)
            fprintf(charge_log,"+%c10",'B'+j);
            break;
         }

      nbp = &neighbour_array[i];      /* for all alpha neighbours */
      for (j=0; j<nbp->n_ligands; j++)
      {
         aap = &mp->atom_array[nbp->atoms[j]];  /* alpha atom */
         abp = &mp->bond_array[nbp->bonds[j]];  /* alpha bond */
                                           /* fetch alpha path conductivity */
         for (k=0, ppa=alpha_path_table; k<nalphapath; k++, ppa++)
              if (AtomSymbolMatch(ap->atom_symbol,
                                   ppa->path.atom_symbol)                &&
                   LigandMatches(aap, abp, &ppa->path.ligands[0], FALSE) &&
                   TRUE)
             break;
         if (k == nalphapath)
         {
            snprintf(msg_buffer, MAXMSG, "%10s: no alpha path for atom %d", mp->name,i+1);
            AddMsgToList(msg_buffer);
            result = FALSE;
            break;
         }

         if (abp->bond_type != SINGLE)
         {
            if (charge_log && old_values[i] != 0)
               fprintf(charge_log,"+%d*%c6",
                                  (abp->bond_type-SINGLE),
                                  'B'+(cip-charge_inc_table));
            ap->value += (abp->bond_type-SINGLE)*cip->mult_inc;
            if (charge_log && old_values[i] != 0) cip->mult_inc_used++;
         }

                                         /* fetch alpha charge increment */
         if (aap->charge != 0)
         {
            for (k=0, cip=charge_inc_table;
                 k<ncharge;
                 k++, cip++)
            if (AtomSymbolMatch(aap->atom_symbol, cip->atom_symbol))
               break;
            if (k == ncharge)
            {
               snprintf(msg_buffer, MAXMSG, "%10s: no alpha increment for atom %d", mp->name,i+1);
               AddMsgToList(msg_buffer);
               result = FALSE;
               break;
            }
            if (charge_log && old_values[i] != 0)
                fprintf(charge_log,"+%d*%c4", aap->charge,'B'+(cip-charge_inc_table));
            ap->value += aap->charge*cip->alpha_inc;
            if (charge_log && old_values[i] != 0) cip->alpha_inc_used++;
         }

                                         /* fetch alpha acidity increment */
         for (k=0, ip=atom_acidity_table; k<natomacidity; k++, ip++)
            if (AtomSymbolMatch(aap->atom_symbol, ip->atom_symbol))
               break;
         if (k == natomacidity)
         {
            snprintf(msg_buffer, MAXMSG, "%10s: no alpha increment for atom %d", mp->name,i+1);
            AddMsgToList(msg_buffer);
            result = FALSE;
            break;
         }
         if (charge_log && old_values[i] != 0)
            fprintf(charge_log,"+%c16", 'B'+(ppa-alpha_path_table));
         if (charge_log && old_values[i] != 0)
            fprintf(charge_log,"*%c11", 'B'+(ip-atom_acidity_table));
               ap->value += ppa->cond*ip->alpha_inc;
         if (charge_log && old_values[i] != 0) ppa->cond_used++;
         if (charge_log && old_values[i] != 0) ip->alpha_inc_used++;

                                    /* for all beta neighbours */
         nbph = &neighbour_array[nbp->atoms[j]];
         for (k=0; k<nbph->n_ligands; k++)
         {
            if (nbph->atoms[k] == i) continue;  /* no loop back */

            bap = &mp->atom_array[nbph->atoms[k]]; /* beta atom */
            bbp = &mp->bond_array[nbph->bonds[k]]; /* beta bond */
            for (l=0, ppb=beta_path_table; l<nbetapath; l++, ppb++) /* fetch beta conductivity */
               if (AtomSymbolMatch(ap->atom_symbol, ppb->path.atom_symbol) &&
                   LigandMatches(aap, abp, &ppb->path.ligands[0], FALSE) &&
                   LigandMatches(bap, bbp, &ppb->path.ligands[1], FALSE) &&
                   TRUE)
                  break;
            if (l == nbetapath)
            {
               snprintf(msg_buffer, MAXMSG, "%10s: no beta increment for atom %d", mp->name,i+1);
               AddMsgToList(msg_buffer);
               result = FALSE;
               break;
            }

                                      /* fetch beta acidity increment */
            for (l=0, ip=atom_acidity_table; l<natomacidity; l++, ip++)
               if (AtomSymbolMatch(bap->atom_symbol, ip->atom_symbol))
               break;
            if (l == natomacidity)
            {
               snprintf(msg_buffer, MAXMSG, "%10s: no beta increment for atom %d", mp->name,i+1);
               AddMsgToList(msg_buffer);
               result = FALSE;
               break;
            }
            if (charge_log && old_values[i] != 0)
               fprintf(charge_log,"+%c20", (int)('B'+(ppb-beta_path_table)));
//   fprintf(charge_log,"+%c16*%c20", (int)('B'+(ppa-alpha_path_table)),
//                                    (int)('B'+(ppb-beta_path_table)));
            if (charge_log && old_values[i] != 0)
               fprintf(charge_log,"*%c12", 'B'+(ip-atom_acidity_table));
//         ap->value += ppa->cond*ppb->cond*ip->beta_inc;
            ap->value += ppb->cond*ip->beta_inc;
            if (charge_log && old_values[i] != 0) ppb->cond_used++;
            if (charge_log && old_values[i] != 0) ip->beta_inc_used++;

                                          /* fetch beta charge increment */
            if (bap->charge != 0)
            {
               for (l=0, cip=charge_inc_table; l<ncharge; l++, cip++)
               if (AtomSymbolMatch(bap->atom_symbol, cip->atom_symbol))
                  break;
               if (l == ncharge)
               {
                  snprintf(msg_buffer, MAXMSG, "%10s: no beta increment for atom %d", mp->name,i+1);
                  AddMsgToList(msg_buffer);
                  result = FALSE;
                  break;
               }
               if (charge_log && old_values[i] != 0)
                  fprintf(charge_log,"+%d*%c5", bap->charge,'B'+(cip-charge_inc_table));
               ap->value += bap->charge*cip->beta_inc;
               if (charge_log && old_values[i] != 0) cip->beta_inc_used++;
            }
         }
      }
      if (charge_log && old_values[i] != 0.0) fprintf(charge_log,"\t%g\n",ap->value);
      ap->value = 7 + beta*(ap->value-7)
                    + alpha*(ap->value-7)*
                    alpha*(ap->value-7)*
                    alpha*(ap->value-7);
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
   int i, result=0;
   struct reaccs_atom_t *ap;
   double epsilon=0.000001;

   for (i=0, min_pKa=1000.0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color != NONE  &&  ap->value < min_pKa)
         min_pKa = ap->value;

   next_pKa = 1000.0;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color != NONE  &&  ap->value < min_pKa+epsilon)
      {
         result++;
         ap->color = NONE+1;
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

   for (i=0; i<mp->n_atoms; i++)
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
   double epsilon=0.000001;

                                             /* set electronegativities */
   atom_ranks = TypeAlloc(mp->n_atoms, struct atom_rank_t);
   for (i=0; i<mp->n_atoms; i++)
   {
      atom_ranks[i].atom = i+1;
      atom_ranks[i].rank = 0;
      atom_ranks[i].rank_sum = 0;
      for (j=0; j<nelneg; j++)
         if (0 == strcmp(mp->atom_array[i].atom_symbol,elneg_table[j].symbol))
         {
            atom_ranks[i].elneg = elneg_table[j].value;
            break;
         }
      if (j == nelneg)
      {
         fprintf(stderr,"atom symbol '%s' not in periodic table\n", mp->atom_array[i].atom_symbol);
         MyFree((char *)atom_ranks);
         return (-1);
      }
      atom_ranks[i].elneg += 3.0*mp->atom_array[i].charge -
                              elneg_table[0].value;
   }

                        /* do preliminary ranking based on el. neg. */
   for (i=1; i<mp->n_atoms; i++)    /* sort by decreasing elneg. */
      for (j=i; j>0; j--)
         if (atom_ranks[j].elneg > atom_ranks[j-1].elneg+epsilon)
         {
            ar_tmp = atom_ranks[j];
            atom_ranks[j] = atom_ranks[j-1];
            atom_ranks[j-1] = ar_tmp;
         }
         else
            break;
   atom_ranks[0].rank = 0;            /* set ranks */
   for (i=1, j=0; i<mp->n_atoms; i++)
   {
      if (atom_ranks[i].elneg < atom_ranks[i-1].elneg-epsilon) j = i;
      atom_ranks[i].rank = j;
   }
   for (i=1; i<mp->n_atoms; i++)  /* resort by atom number */
      for (j=i; j>0; j--)
         if (atom_ranks[j].atom < atom_ranks[j-1].atom)
         {
            ar_tmp = atom_ranks[j];
            atom_ranks[j] = atom_ranks[j-1];
            atom_ranks[j-1] = ar_tmp;
         }
         else
            break;

       /* use unsaturation to split ranks (rank sum is misused here) */
   for (i=0; i<mp->n_atoms; i++) atom_ranks[i].rank_sum = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->bond_type != SINGLE)
      {
         atom_ranks[bp->atoms[0]-1].rank_sum++;
         atom_ranks[bp->atoms[1]-1].rank_sum++;
      }
   for (i=1; i<mp->n_atoms; i++) /* sort by rank(asc.) + unsat.(desc.) */
      for (j=i; j>0; j--)
         if (atom_ranks[j].rank < atom_ranks[j-1].rank  ||
           (atom_ranks[j].rank == atom_ranks[j-1].rank  &&
            atom_ranks[j].rank_sum > atom_ranks[j-1].rank_sum))
         {
            ar_tmp = atom_ranks[j];
            atom_ranks[j] = atom_ranks[j-1];
            atom_ranks[j-1] = ar_tmp;
         }
         else
            break;
   for (i=1, j=0; i<mp->n_atoms; i++) /* set new ranks */
   {
      if (atom_ranks[i].rank > atom_ranks[i-1].rank  ||
          atom_ranks[i].rank_sum < atom_ranks[i-1].rank_sum) j = i;
         atom_ranks[i].rank = j;
   }
   for (i=1; i<mp->n_atoms; i++)    /* restore atom number order */
      for (j=i; j>0; j--)
         if (atom_ranks[j].atom < atom_ranks[j-1].atom)
         {
            ar_tmp = atom_ranks[j];
            atom_ranks[j] = atom_ranks[j-1];
            atom_ranks[j-1] = ar_tmp;
       }
       else
          break;

   for (i=0; i<mp->n_atoms; i++) atom_ranks[i].n_ligands = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      atom_ranks[bp->atoms[0]-1].n_ligands++;
      atom_ranks[bp->atoms[1]-1].n_ligands++;
   }

                         /* refine ranking using neighbour rank sums */
   do_cis_trans = FALSE;
   do
   {
      for (i=0; i<mp->n_atoms; i++)
         numbering[i] = atom_ranks[i].rank;
                             /* compute rank sums */
      for (i=0; i<mp->n_atoms; i++) atom_ranks[i].rank_sum = 0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         atom_ranks[bp->atoms[0]-1].rank_sum +=
         atom_ranks[bp->atoms[1]-1].rank;
         atom_ranks[bp->atoms[1]-1].rank_sum +=
         atom_ranks[bp->atoms[0]-1].rank;
      }
      if (do_cis_trans)
      {
         CisTransPerception(mp, numbering);
         for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
            if (bp->color != NONE)
            {
               atom_ranks[bp->atoms[0]-1].rank_sum += bp->color;
               atom_ranks[bp->atoms[1]-1].rank_sum += bp->color;
            }
      }
      for (i=0; i<mp->n_atoms; i++) /* use average rank sum */
         if (atom_ranks[i].n_ligands > 0)
         {
            atom_ranks[i].rank_sum *= 10; /* shift dec. point */
            atom_ranks[i].rank_sum /= atom_ranks[i].n_ligands;
         }
      for (i=1; i<mp->n_atoms; i++)   /* sort by rank + ranksum */
         for (j=i; j>0; j--)
            if (atom_ranks[j].rank < atom_ranks[j-1].rank  ||
                (atom_ranks[j].rank == atom_ranks[j-1].rank  &&
                 atom_ranks[j].rank_sum < atom_ranks[j-1].rank_sum))
            {
               ar_tmp = atom_ranks[j];
               atom_ranks[j] = atom_ranks[j-1];
               atom_ranks[j-1] = ar_tmp;
            }
            else
               break;
                           /* compute new ranks */
      for (i=1, changed = FALSE, j=0; i<mp->n_atoms; i++)
      {
         if (atom_ranks[i].rank > atom_ranks[i-1].rank  ||
             atom_ranks[i].rank_sum > atom_ranks[i-1].rank_sum) j = i;
         changed = changed || atom_ranks[i].rank != j;
         atom_ranks[i].rank = j;
      }
      for (i=1; i<mp->n_atoms; i++) /* restore atom number order */
         for (j=i; j>0; j--)
            if (atom_ranks[j].atom < atom_ranks[j-1].atom)
            {
               ar_tmp = atom_ranks[j];
               atom_ranks[j] = atom_ranks[j-1];
               atom_ranks[j-1] = ar_tmp;
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
   for (i=0, min_rank=mp->n_atoms; i<mp->n_atoms; i++)
      if (mp->atom_array[i].color != NONE  &&
          atom_ranks[i].rank < min_rank)
         min_rank = atom_ranks[i].rank;
   for (i=0; i<mp->n_atoms; i++)
   {
      if (mp->atom_array[i].color != NONE  &&
          atom_ranks[i].rank == min_rank)
      {                        /* count members of minimum class */
         result++;
         mp->atom_array[i].value = atom_ranks[i].rank+1;
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
      for (i=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].color != NONE)
         {
            sprintf(msg_buffer,"atom %d in minimal rank class",i+1);
            AddMsgToList(msg_buffer);
         }
   }
   return (result);
}

static double acidity_limit = 24.0;

void SetAcidityLimit(double limit)
/*
 * Sets the acidity limit for atoms considered to be available
 * for deprotonation.
 */
{
   acidity_limit = limit;
}

int RechargeMolecule(struct reaccs_molecule_t *mp,
                     int desired_charge,
                     int *ndeprot, int *nrefine)
/*
 * Removes hydrogens from *mp until desired_charge is reached. The
 * positions for hydrogen removal are selected by "acidity" combined
 * with a refinement algorithm. It returns TRUE if molecule could be
 * neutralized and FALSE if any problem were encountered.
 * *ndeprot and *nrefine are set to the number of deprotonations
 * and refinement cycles performed.
 */
{
   int i;
   int nacid;
   double gap, pKa_value;
   int *numbering;

   *ndeprot = 0;   /* number of deprotonation cycles */
   *nrefine = 0;    /* number of refinements necessary */
   for (;;)
   {
      while (TotalCharge(mp) > desired_charge)
      {
         SetpKaValues(mp);
         nacid = MarkMostAcidicAtoms(mp, &pKa_value, &gap);
         if (pKa_value > acidity_limit)
         {
            sprintf(msg_buffer,"pKa_value = %.2g",pKa_value);
            AddMsgToList(msg_buffer);
            return (FALSE);     /* not acidic enough */
         }
         else if (nacid == 1 ||               /* easy case */
                  (nacid == TotalCharge(mp) && gap > 8.0))
         {
            sprintf(msg_buffer,"pKa = %.2g",pKa_value);
            AddMsgToList(msg_buffer);
            DecrementMarkedCharges(mp);
         }
         else
            break;
         (*ndeprot)++;
      }

      if (TotalCharge(mp) > desired_charge)
      {
         (*nrefine)++;
         numbering = TypeAlloc(mp->n_atoms, int);
         nacid = RefineAcidicAtoms(mp, numbering);
         if (nacid == 1)
            DecrementMarkedCharges(mp);
         else if (AllCentersRefined(mp, numbering))
         {
            for (i=0; i<mp->n_atoms; i++)       /* Select one mark */
               if (mp->atom_array[i].color != NONE && nacid != 1)
               {
                  nacid--; mp->atom_array[i].color = NONE;
               }
            DecrementMarkedCharges(mp);
         }
         else
         {
            snprintf(msg_buffer, MAXMSG, "%10s: could not fix charges", mp->name);
            AddMsgToList(msg_buffer);
            return (FALSE);
         }
         MyFree((char *)numbering);
         (*ndeprot)++;
      }
      else
      {
         ResetValues(mp);
         return (TRUE);
      }
   }
}
