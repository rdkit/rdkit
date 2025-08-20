/**********************************************************************/
/*                                                                    */
/*    File:           matchtest.c                                     */
/*                                                                    */
/*    Purpose:        This file implements the functions needed to    */
/*                    test the substructure matching components.      */
/*                                                                    */
/*    History:        11-Oct-2012     Start of development.           */
/*                                                                    */
/**********************************************************************/

#ifdef __WIN32__
#include <windows.h>
#endif

#include "local.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "reaccs.h"
#include "utilities.h"
#include "perceive.h"
#include "reaccsio.h"
#include "stereo.h"
#include "pattern.h"

#include "smi2mol.h"
#include "canonizer.h"

#include "ssmatch.h"

struct generator_t
{
   int use_generator;
   char *name;
   char *short_name;
   unsigned int mask;
   int *counts;
}
generators[] =
{
   {FALSE, "RING_PATTERN", "RPT", USE_RING_PATTERN, (int *)NULL},
   {FALSE, "RING_PATH", "RPH", USE_RING_PATH, (int *)NULL},
   {FALSE, "ATOM_SYMBOL_PATH", "ASP", USE_ATOM_SYMBOL_PATH, (int *)NULL},
   {FALSE, "ATOM_CLASS_PATH", "ACP", USE_ATOM_CLASS_PATH, (int *)NULL},
   {FALSE, "ATOM_COUNT", "AC", USE_ATOM_COUNT, (int *)NULL},
   {FALSE, "AUGMENTED_ATOM", "AA", USE_AUGMENTED_ATOM, (int *)NULL},
   {FALSE, "HCOUNT_PATH", "HPH", USE_HCOUNT_PATH, (int *)NULL},
   {FALSE, "HCOUNT_CLASS_PATH", "HCP", USE_HCOUNT_CLASS_PATH, (int *)NULL},
   {FALSE, "HCOUNT_PAIR", "HP", USE_HCOUNT_PAIR, (int *)NULL},
   {FALSE, "BOND_PATH", "BP", USE_BOND_PATH, (int *)NULL},
   {FALSE, "AUGMENTED_BOND", "AB", USE_AUGMENTED_BOND, (int *)NULL},
   {FALSE, "RING_SIZE_COUNTS", "RSC", USE_RING_SIZE_COUNTS, (int *)NULL},
   {FALSE, "DEGREE_PATH", "DP", USE_DEGREE_PATH, (int *)NULL},
   {FALSE, "CLASS_SPIDERS", "CS", USE_CLASS_SPIDERS, (int *)NULL},
   {FALSE, "FEATURE_PAIRS", "FP", USE_FEATURE_PAIRS, (int *)NULL},
   {FALSE, "SCAFFOLD_IDS", "SI", USE_SCAFFOLD_IDS, (int *)NULL},
   {FALSE, "SCAFFOLD_COLORS", "SC", USE_SCAFFOLD_COLORS, (int *)NULL},
   {FALSE, "SCAFFOLD_LINKS", "SL", USE_SCAFFOLD_LINKS, (int *)NULL},
   {FALSE, "SHORTCUT_LABELS", "SCL", USE_SHORTCUT_LABELS, (int *)NULL},
};
int selected_generators = FALSE;

#define MXLN        2000

int as_query = FALSE;
int add_titles = FALSE;
#define BIN_OUTPUT 0
#define LST_OUTPUT 1
#define CNT_OUTPUT 2
int output_format = BIN_OUTPUT;
int nbits = 512;
int lcnt = -1;
int ucnt = -1;
char* input_file;
char src_title[MXLN];

static struct reaccs_molecule_t *SMIToMOLNoStereo(const char *const_smiles)
{
   char *smiles, *cp;
   struct reaccs_molecule_t *mp;
   struct reaccs_bond_t *bp;
   int *H_count;
   neighbourhood_t *nbp;
   int i;

   if (const_smiles)
   {
      smiles = TypeAlloc(strlen(const_smiles)+1, char);
      cp = smiles;
      do
      {
         if ((*const_smiles) != '@') (*cp++) = (*const_smiles);
         const_smiles++;
      }
      while (*const_smiles);
      mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
      MyFree(smiles);
      if (mp == (struct reaccs_molecule_t *)NULL)
      {
         return (mp);
      }
      nbp = TypeAlloc(mp->n_atoms, neighbourhood_t);
      SetupNeighbourhood(mp,nbp,mp->n_atoms);
      SetRingSizeFlags(mp, 14, nbp);
      MyFree((char *)nbp);
      /* Set up hydrogen count fields in structure for matching */
      H_count = TypeAlloc(mp->n_atoms+1, int);
      ComputeImplicitH(mp, H_count);
// fprintf(stderr, "after ComputeImplicitH\n"); Fortify_ListAllMemory();
      /* Add the explicit hydrogens to the implicit counts */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            H_count[bp->atoms[1]]++;
         else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
            H_count[bp->atoms[0]]++;
      }
      /* set the 'query_H_count' field to the correct value */
      for (i=0; i<mp->n_atoms; i++)
         if (H_count[i+1] >= 0) 
            mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];
      MyFree((char *)H_count);
// fprintf(stderr, "after MyFree(H_count)\n"); Fortify_ListAllMemory();
   }
   else
      mp = (struct reaccs_molecule_t *)NULL;
   return (mp);
}

void CountFPFrequency(char *infile)
{
   int i, j, ncounts;
   FILE *fp;
   char line[MXLN];
   char smiles[2*MXLN];
   int nsmi;
   struct reaccs_molecule_t *mp;
   int* fp_counts;

   fp = fopen(infile, "r");
   if (IsNULL(fp))
   {
      fprintf(stderr, "Could not open '%s'\n", infile);
      exit (EXIT_FAILURE);
   }
   fp_counts = TypeAlloc(nbits, int);
   for (i=0; i<sizeof(generators)/sizeof(struct generator_t); i++)
      if (generators[i].use_generator)
      {
         generators[i].counts = TypeAlloc(nbits, int);
      }
   nsmi = 0;
   while (!feof(fp))
   {

      fgets(line, MXLN-1, fp);
      smiles[0] = '\0';
      sscanf(line, "%s", smiles);
      mp = SMIToMOLNoStereo(smiles);
      if (nsmi == 0  &&  IsNULL(mp))
      {
         strncpy(src_title, line, MXLN-1); src_title[MXLN-1] = '\0';
         continue;
      }
      // here we have a regular molecule
      nsmi++;
      for (i=0; i<sizeof(generators)/sizeof(struct generator_t); i++)
      {
         if (generators[i].use_generator)
         {
            SetFingerprintCountsWithFocus(mp, fp_counts, nbits, generators[i].mask, as_query, USE_DY_AROMATICITY, 0);
            for (j=0; j<nbits; j++)
               if (fp_counts[j] > 0)
                  generators[i].counts[j]++;
         }
      }
      FreeMolecule(mp);
      if (feof(fp)) break;

   }
   ncounts = 0;
   for (i=0; i<sizeof(generators)/sizeof(struct generator_t); i++)
      for (j=0; j<nbits; j++)
         if (generators[i].use_generator  &&  generators[i].counts[j] > 0) ncounts++;
   MyFree((char *)fp_counts);
   fprintf(stderr, "%d molecules counted. Found %d bits\n", nsmi, ncounts);
   fclose(fp);
}

static void explain(char * prognam, int argc)
{
    if (strstr(prognam, "matchtest"))
    {
       fprintf(stderr,"matchtest used with %d arguments\n",argc-1);
       fprintf(stderr,"usage: matchtest <input-smiles-file> {<query-mol-file>}\n");
       fprintf(stderr,"Reads the SMILES in the first file argument and test SS matches to the subsequent MOL file arguments\n");
       exit (EXIT_FAILURE);
    }
    else if (strstr(prognam, "mol2fp"))
    {
       fprintf(stderr,"usage: mol2fp {options} <mol-file>\n");
       fprintf(stderr,"Reads the gigen MOL file and returns fingerprints according to the options\n");
       fprintf(stderr,"The following options are defined:\n");
       fprintf(stderr,"-q\t                \tFingerprint MOL file as query\n");
       fprintf(stderr,"-t\t                \tAdd column titles\n");
       fprintf(stderr,"-f\tbin|lst|cnt     \toutput format (0/1, list, or feature counts)\n");
       fprintf(stderr,"-g\t<generator-list>\tComma separated list of FP generators, pick from\n");
       fprintf(stderr,"  \t                \tRING_PATTERN, RING_PATH, ATOM_SYMBOL_PATH, ATOM_CLASS_PATH, ATOM_COUNT,\n");
       fprintf(stderr,"  \t                \tAUGMENTED_ATOM, HCOUNT_PATH, HCOUNT_CLASS_PATH, BOND_BATH, AUGMENTED_BOND,\n");
       fprintf(stderr,"  \t                \tRING_SIZE_COUNTS, DEGREE_PATH, CLASS_SPIDERS, FEATURE_PAIRS,\n");
       fprintf(stderr,"  \t                \t(non-SS-features: SCAFFOLD_COLORS, SCAFFOLD_LINKS, SHORTCUT_LABELS not used for -q)\n");
       fprintf(stderr,"-n\tnbits           \tWidth of fingerprint in bits (will be rounded up to multiples of 8)\n");
       fprintf(stderr,"  \t                \tMissing -g option causes all generators to fire anonymously\n");
       exit (EXIT_FAILURE);
    }
    else if (strstr(prognam, "smi2fp"))
    {
       fprintf(stderr,"usage: smi2fp {options} <smiles-input-file> ><fp-augmented-output-file>\n");
       fprintf(stderr,"The following options are defined:\n");
       fprintf(stderr,"-q\t                \tFingerprint SMILES strings as queries\n");
       fprintf(stderr,"-t\t                \tAdd column titles\n");
       fprintf(stderr,"-f\tbin|lst|cnt     \toutput format (0/1, list, or feature counts)\n");
       fprintf(stderr,"-g\t<generator-list>\tComma separated list of FP generators, pick from\n");
       fprintf(stderr,"  \t                \tRING_PATTERN, RING_PATH, ATOM_SYMBOL_PATH, ATOM_CLASS_PATH, ATOM_COUNT,\n");
       fprintf(stderr,"  \t                \tAUGMENTED_ATOM, HCOUNT_PATH, HCOUNT_CLASS_PATH, BOND_BATH, AUGMENTED_BOND,\n");
       fprintf(stderr,"  \t                \tRING_SIZE_COUNTS, DEGREE_PATH, CLASS_SPIDERS, FEATURE_PAIRS,\n");
       fprintf(stderr,"  \t                \t(non-SS-features: SCAFFOLD_COLORS, SCAFFOLD_LINKS, SHORTCUT_LABELS not used for -q)\n");
       fprintf(stderr,"-n\tnbits           \tWidth of fingerprint in bits (will be rounded up to multiples of 8)\n");
       fprintf(stderr,"-l\t<min-count>     \tIgnore FP bits with less than <min-count> occurences in input file\n");
       fprintf(stderr,"-u\t<max-count>     \tIgnore FP bits with more than <max-count> occurences in input file\n");
       fprintf(stderr,"  \t                \t-l and -u cause a pre-scan of the input to collect the sample counts\n");
       fprintf(stderr,"  \t                \tMissing -g option causes all generators to fire anonymously\n");
       exit (EXIT_FAILURE);
    }
}

#define MXBUF	20000

static void init_options(int argc, char *argv[])
{
   int i, j;
   int itmp;
   char buffer[MXBUF];
   for (i=1; i<argc; i++)
   {
      if (0 == strcmp("-q", argv[i]))
      {
         as_query = TRUE;
         fprintf(stderr, "fingerprinting as query\n");
         continue;
      }
      if (0 == strcmp("-t", argv[i]))
      {
         add_titles = TRUE;
         fprintf(stderr, "adding titles\n");
         continue;
      }
      if (0 == strcmp("-f", argv[i]))
      {
         i++;
         if (0 == strcmp("bin", argv[i]))
         {
            output_format = BIN_OUTPUT;
            fprintf(stderr, "output format BIN\n");
         }
         else if (0 == strcmp("lst", argv[i]))
         {
            output_format = LST_OUTPUT;
            fprintf(stderr, "output format LST\n");
         }
         else if (0 == strcmp("cnt", argv[i]))
         {
            output_format = CNT_OUTPUT;
            fprintf(stderr, "output format CNT\n");
         }
         else
         {
            fprintf(stderr, "Error setting output format -f\n");
            explain(argv[0], argc);
         }
         continue;
      }
      if (0 == strcmp("-g", argv[i]))
      {
         selected_generators = TRUE;
         i++;
         for (j=0; j<sizeof(generators)/sizeof(struct generator_t); j++)
            if (!IsNULL(strstr(argv[i], generators[j].name)))
            {
               generators[j].use_generator = TRUE;
               fprintf(stderr, "found genrator '%s' in list %s\n", generators[j].name, argv[i]);
            }
         continue;
      }
      if (0 == strcmp("-n", argv[i]))
      {
         i++;
         itmp = atoi(argv[i]);
         for (nbits=1; nbits<itmp; nbits *= 2)
            ;
         fprintf(stderr, "fingerprint width %d bits computed from '-n %d'\n", nbits, itmp);
         continue;
      }
      if (0 == strcmp("-l", argv[i]))
      {
         i++;
         lcnt = atoi(argv[i]);
         fprintf(stderr, "fingerprint lower count limit = %d\n", lcnt);
         continue;
      }
      if (0 == strcmp("-u", argv[i]))
      {
         i++;
         ucnt = atoi(argv[i]);
         fprintf(stderr, "fingerprint upper count limit = %d\n", ucnt);
         continue;
      }
      if (i+1 < argc)
      {
         fprintf(stderr, "Illegal option\n");
         explain(argv[0], argc);
      }
      input_file = argv[i];
      fprintf(stderr, "reading from file = %s\n", input_file);
      return;    // no further options => grab input and output files
   }
}

#define MAXLINE 20000

#include "didepict.h"

/* used for 16-bit Borland Turbo C compilation */
#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

FILE *log_file;

static struct reaccs_molecule_t *QueryFromFile(char *fname)
{
   struct reaccs_molecule_t *qp;
   Fortran_FILE *ffp;
   neighbourhood_t *nbp;

   ffp = FortranOpen(fname,"r");
   if (IsNULL(ffp)) return (struct reaccs_molecule_t *)NULL;
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read query file '%s'\n", fname);
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);

   if (qp == (struct reaccs_molecule_t *)NULL) return qp;
   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   qp = ThreadQueryAtoms(qp);
   return (qp);
}

extern int refine_success;

int main(int argc, char *argv[])
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   char *smiles;
   char line[MAXLINE];
   char buffer[MXBUF];
   char *cp, *smistart;
   char prognam[256];
   int nsmi, iquery, nmatch, nfail, nrefine;
   FILE *infile;
   ssmatch_t * matches;

   int i, j, jj, ncols, total_cols;
   int* fp_counts, *fp_counts_tmp;
   int result;

// FORTIFY long bytes_allocated;
// FORTIFY Fortify_SetOutputFunc(logFortifyMessage);
// FORTIFY Fortify_OutputStatistics();
// FORTIFY fprintf(stderr, "fortify enabled\n");

   log_file = stderr;

   strncpy(prognam, argv[0], 255);
   for (cp=prognam; (*cp); cp++)
      (*cp) = tolower(*cp);

   /* Stand-alone driver for substructure matching */
   if (strstr(prognam, "matchtest"))
   {
      if( argc < 2 )
      {
         explain(argv[0], argc);
         return (EXIT_FAILURE);
      }
// FORTIFY Fortify_OutputStatistics();
      SetStripZeroes(TRUE);
      for (iquery=2; iquery<argc; iquery++)
      {
         nsmi = 0;
         nmatch = 0; nfail = 0; nrefine = 0;
         // read query MOL file
         qp = QueryFromFile(argv[iquery]);
         if (qp == (struct reaccs_molecule_t *)NULL) continue;
         // step through all SMILES
         infile = fopen(argv[1], "r");
         while (!feof(infile))
         {
            fgets(line, MAXLINE-1, infile);
            if (strstr(line, "EXIT")) break;
            /* remove any remaining '\n' */
            line[MAXLINE-1] = '\0';
            if (strlen(line) > 0 && line[strlen(line)-1] == '\n')
               line[strlen(line)-1] = '\0';
            smistart = line;
            /* get SMILES as first word in line */
            for (cp=smistart; (*cp) && !isspace(*cp); cp++)
               ;
            strncpy(buffer, smistart, cp-smistart);
            buffer[cp-smistart] = '\0';
            if (feof(infile)) break;
            if (strlen(buffer) == 0) continue;
            mp = SMIToMOLNoStereo(buffer);
            if (mp == (struct reaccs_molecule_t *)NULL) continue;
// FORTIFY Fortify_EnterScope();
// FORTIFY bytes_allocated = fortifySet();
            // convert SMILES to molecule
// fprintf(stdout, "%s\n", buffer);
            // PrintREACCSMolecule(stderr, mp, "");
            // do the match
// fprintf(stderr, "before SSMatch\n"); Fortify_ListAllMemory();
            matches = SSMatch(mp, qp, TRUE);
// fprintf(stderr, "after SSMatch\n"); Fortify_ListAllMemory();
            if (matches != (ssmatch_t *)NULL)
            {
               nmatch++;
               FreeSSMatch(matches);
               FreeSSMatchHeap();
// fprintf(stderr, "after FreeSSMatchHeap\n"); Fortify_ListAllMemory();
            }
            else
               nfail++;
            if (refine_success)
            {
               nrefine++;
// fprintf(stdout, "%s\n", buffer);
            }
            // free memory
// FORTIFY fortifyTest(bytes_allocated, buffer);
// FORTIFY Fortify_LeaveScope();
            FreeMolecule(mp);
         }
         fclose(infile);
         // return query memory
         FreeMolecule(qp);
         fprintf(stderr, "%s:\t%d matches\t%d non-matches\t%d refine-matches\n", argv[iquery], nmatch, nfail, nrefine);
      }
// FORTIFY Fortify_OutputStatistics();
// FORTIFY fprintf(stderr, "fortify checked\n");
      return (EXIT_SUCCESS);
   }

   // collect the fingerprinting parameters from the command line for both 'mol2fp' and 'smi2fp'

   init_options(argc, argv);

   // perform the fingerprint computations
   if (strstr(prognam, "mol2fp"))  // print annotated fingerprints of size 65536
   {
      if(IsNULL(infile))
      {
         explain(argv[0], argc);
         return (EXIT_FAILURE);
      }
      // read query MOL file
      qp = QueryFromFile(input_file);
      if (qp == (struct reaccs_molecule_t *)NULL)
      {
         fprintf(stderr, "Could not read query structure %s\n", argv[1]);
         return (EXIT_FAILURE);
      }
      fprintf(stderr, "Processing query '%s'\n", argv[1]);
      fp_counts = TypeAlloc(8192*8, int);
      fp_counts_tmp = TypeAlloc(8192*8, int);
      for (i=0; i<32; i++)
      {
         if (IsNULL(fp_prefixes[i])) break;
         result = SetFingerprintCountsWithFocus(qp, fp_counts, 8192*8, (1<<i), TRUE, USE_DY_AROMATICITY, 0);
         if (result <= 0) continue;
         fprintf(stderr, "%s:\t", fp_prefixes[i]);
         for (j=0; j<8192*8; j++)
            if (fp_counts[j] > 0)
            {
               fprintf(stderr, "%d@",j);
               for (jj=0; jj<qp->n_atoms; jj++)
               {
                  result = SetFingerprintCountsWithFocus(qp, fp_counts_tmp, 8192*8, (1<<i), TRUE, USE_DY_AROMATICITY, jj+1);
                  if (fp_counts[j]!=fp_counts_tmp[j])
                  {
                     fprintf(stderr, "%s(%d)", qp->atom_array[jj].atom_symbol, jj+1);
                  }
               }
               fprintf(stderr, "\t");
            }
         fprintf(stderr, "\n");
      }
   }

   else if (strstr(prognam, "smi2fp"))  // print annotated fingerprints of size 65536
   {
      if (selected_generators) CountFPFrequency(input_file);
      if (0 == strcmp("-",input_file))
         infile = stdin;
      else
      {
         infile = fopen(input_file, "r");
         if (IsNULL(infile))
         {
            fprintf(stderr, "Could not open '%s'\n", input_file);
            exit (EXIT_FAILURE);
         }
      }
      if (add_titles)    // extend the first line of the data file and extend it by the fingerprint titles
      {
         fprintf(stderr, "Adding title, lcnt=%d, ucnt=%d\n", lcnt, ucnt);
         fgets(line, MAXLINE-1, infile);
         line[strlen(line) - 1] = '\0';
         fprintf(stdout, "%s", line);
         if (selected_generators)
         {
            total_cols = 0;
            for (i=0; i<sizeof(generators)/sizeof(struct generator_t); i++)
               if (generators[i].use_generator)
               {
                  ncols = 0;
                  for (j=0; j<nbits; j++)
                     if (generators[i].counts[j] > lcnt  &&  (ucnt < 0  ||  ucnt > generators[i].counts[j]))
                     {
                        fprintf(stdout, "\t%s(%d)", generators[i].short_name, j);
                        ncols++;
                     }
                  fprintf(stderr, " %d FP columns for generator '%s'\n", ncols, generators[i].name);
                  total_cols += ncols;
               }
            fprintf(stderr, " %d FP columns for all generators\n", total_cols);
         }
         else
         {
         }
         fprintf(stdout, "\n");
      }
      nsmi = 0;
      fp_counts = TypeAlloc(nbits, int);
      while (!feof(infile))
      {
         if (++nsmi%10 == 0) fprintf(stderr,".");
         if (nsmi%500 == 0) fprintf(stderr," %d\n", nsmi);
         if (nsmi > 10000) exit(EXIT_FAILURE);

         fgets(line, MXLN-1, infile);
         line[strlen(line) - 1] = '\0';
         buffer[0] = '\0';
         sscanf(line, "%s", buffer);
         mp = SMIToMOLNoStereo(buffer);
         if (IsNULL(mp)) continue;
         // here we have a regular molecule
         fprintf(stdout, "%s", line);
         for (i=0; i<sizeof(generators)/sizeof(struct generator_t); i++)
         {
            if (generators[i].use_generator)
            {
               SetFingerprintCountsWithFocus(mp, fp_counts, nbits, generators[i].mask, as_query, USE_DY_AROMATICITY, 0);
               for (j=0; j<nbits; j++)
                  if (generators[i].counts[j] > lcnt  &&  (ucnt < 0  ||  ucnt > generators[i].counts[j]))
                     fprintf(stdout, "\t%d", fp_counts[j]);
            }
         }
         fprintf(stdout, "\n");
         FreeMolecule(mp);
         if (feof(infile)) break;
      }
      if (nsmi%500 != 0) fprintf(stderr," %d\n", nsmi);
      MyFree((char *)fp_counts);
   }

   else
   {
      explain(prognam, argc);
   }

   return (EXIT_SUCCESS);
}
