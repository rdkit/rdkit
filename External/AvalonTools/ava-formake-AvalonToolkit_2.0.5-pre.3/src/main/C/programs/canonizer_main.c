/**********************************************************************/
/*                                                                    */
/*    File:           canonizer_main.c                                */
/*                                                                    */
/*    Purpose:        This file implements main execeutable to run    */
/*                    the canonicalizer to generate canonical SMILES  */
/*                    from an MDL MOL data structure.                 */
/*                                                                    */
/*    History:        05-Feb-2020     Extracted from canonizer.c.     */
/*                                                                    */
/**********************************************************************/

#ifdef __WIN32__
#include <windows.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "canonizer.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "smi2mol.h"
#include "utilities.h"

static int CheckHeap(void)
{
   int result = TRUE;

/* non-dummy only for 32-bit Windows */
#ifdef __WIN32__
   HGLOBAL heap;

   heap = GetProcessHeap();
   result = HeapValidate(heap,0,NULL);
   if (!result)
      fprintf(stderr, "HeapValidate(%d, %d, %d) returns %d\n", heap, 0, NULL, result);
#endif

   return (result);
}

static void explain(char * prognam, int argc)
{
    if (strstr(prognam, "canonizer"))
   {
      fprintf(stderr,"canonizer used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: canonizer {flags} <input.tbl >output.tbl\n");
      fprintf(stderr,"The following flags are known: \n");
      fprintf(stderr,"-h      \tShows this list \n");
      fprintf(stderr,"-hasid  \tfirst token is an id and appended to output\n");
      fprintf(stderr,"-add    \tappends the result to the input line\n");
      fprintf(stderr,"-unique \talso outputs the SMILES w/o stereo\n");
      fprintf(stderr,"-monomer\tconvert the SMILES to a monomer definition\n");
      fprintf(stderr,"-keep_mapping\tdoes not clear away atom-mapping information\n");
      fprintf(stderr,"-mwmf   \tcalculate MW and MF\n");
      fprintf(stderr,"-d \tuse Daylight aromaticity perception\n");
      fprintf(stderr,"-do_test \tperform a regression test \n");
      exit (EXIT_FAILURE);
   }
}

#include "didepict.h"

/* used for 16-bit Borland Turbo C compilation */
#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

#define MXBUF	MAXLINE

FILE *log_file;

int main(int argc, char *argv[])
{
   struct reaccs_molecule_t *mp;
   char *smiles, *usmiles;
   char *coordp;
   char line[MAXLINE];
   char buffer[MXBUF];
   char check_smiles[MXBUF];
   char *cp;
   char prognam[256];
   int nsmi, i;
   double mw;
   char mf[MAXDATA+1];

   int show_help    = FALSE;
   int add_result   = FALSE;
   int unique_smiles= FALSE;
   int first_run, nruns;
   int do_test      = FALSE;
   int mwmf         = FALSE;
   int hasid        = FALSE;
   char *smistart;
   char id[32];
   int flags        = 0;

   log_file = stderr;

   flags |= REMOVE_MAPPING;
   /* Read options */
   while (argc > 1  &&  argv[1][0] == '-' && argv[1][1] != '\0' )
   {
      if (0 == strcmp(argv[1], "-add"))      add_result = TRUE;
      if (0 == strcmp(argv[1], "-unique"))   unique_smiles = TRUE;
      if (0 == strcmp(argv[1], "-do_test"))  do_test   = TRUE;
      if (0 == strcmp(argv[1], "-hasid"))    hasid   = TRUE;
      if (0 == strcmp(argv[1], "-mwmf"))     mwmf   = TRUE;
      if (0 == strcmp(argv[1], "-d"))       flags   |= DY_AROMATICITY;
      if (0 == strcmp(argv[1], "-monomer"))  flags   |= TO_MONOMER;
      if (0 == strcmp(argv[1], "-keep_mapping"))  flags   &= ~REMOVE_MAPPING;
      if (0 == strcmp(argv[1], "-h"))        show_help  = TRUE;
      for (i=2; i<argc; i++) argv[i-1] = argv[i];
      argc--;
   }

   coordp = NULL;
   strncpy(prognam, argv[0], 255);
   for (cp=prognam; (*cp); cp++)
      (*cp) = tolower(*cp);

   if( show_help )
   {
      explain(prognam, argc);
      return (EXIT_FAILURE);
   }

   /* Stand-alone driver for canonicalization */
   if (strstr(prognam, "canonizer"))
   {
      nsmi = 0;
      SetStripZeroes(TRUE);
      while (!feof(stdin))
      {
         fgets(line, MAXLINE-1, stdin);
if (strstr(line, "EXIT")) break;
         /* remove any remaining '\n' */
         line[MAXLINE-1] = '\0';
         if (strlen(line) > 0 && line[strlen(line)-1] == '\n')
            line[strlen(line)-1] = '\0';
	 if (hasid)
	 {
	    for (cp=line; (*cp) && !isspace(*cp); cp++)
	       ;
	    strncpy(id, line, cp-line);
	    id[cp-line] = '\0';
	    while (isspace(*cp)) cp++;
	    smistart = cp;
	 }
	 else
	    smistart = line;
         /* get SMILES as first word in line */
         for (cp=smistart; (*cp) && !isspace(*cp); cp++)
            ;
         strncpy(buffer, smistart, cp-smistart);
         buffer[cp-smistart] = '\0';
         if (feof(stdin)) break;
         if (strlen(buffer) == 0) continue;
         first_run = TRUE;
         nruns = 0;
// if (do_test) fprintf(stderr, "step 1\n");
         // if (do_test) fprintf(stderr, "%s\n", buffer);

         SetSeed(13);

repeat_do_test:
// if (do_test) fprintf(stderr, "step 1a\n");
         if (do_test)
            smiles = CanSmiles(buffer,
                               SCRAMBLE |
                               flags |
                               DB_STEREO |
                               CENTER_STEREO);
         else
            smiles = CanSmiles(buffer,
                               flags |
                               DB_STEREO |
                               CENTER_STEREO);
// if (do_test) fprintf(stderr, "step 2\n");
         if (unique_smiles)
            usmiles = CanSmiles(buffer, flags);
         else
            usmiles=(char *)NULL;
         if (!smiles)
         {
            /* No result => just print line */
            if (first_run  &&  add_result) fprintf(stdout,"%s\n", line);
            continue;
         }
// if (do_test) fprintf(stderr, "step 3\n");
         nruns++;
         if (first_run)
         {
            nsmi++;
            strcpy(check_smiles, smiles);
         }

// if (do_test) fprintf(stderr, "step 4\n");
         if (add_result)
         {
            if (first_run)
            {
               fprintf(stdout, "%s\t%s", smistart, smiles);
               if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
               fprintf(stdout, "\n");
            }
            else
            {
               if (0 != strcmp(check_smiles, smiles))
               {
                  fprintf(stderr, "'%s' <> '%s'\n", check_smiles, smiles);
                  fprintf(stderr, "original line = %s\n", line);
                  // exit(1);
               }
            }
         }
         else
         {
            if (first_run)
            {
               if (mwmf)
               {
                  mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
                  if (!mp) fprintf(stdout,"%s\n", smistart);
                  else
                  {
                     fprintf(stdout,"%s", smistart);
                     mw = MolecularWeight(mp);
                     MolecularFormula(mp, mf, MAXDATA);
                     fprintf(stdout,"\t%g\t%s",mw,mf);
                     FreeMolecule(mp);
                     if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
		     if (hasid)
			fprintf(stdout, "\t%s\n", id);
		     else
			fprintf(stdout, "\n");
                  }
               }
               else
               {
                  fprintf(stdout, "%s", smiles);
                  if (unique_smiles) fprintf(stdout, "\t%s", usmiles);
		  if (hasid)
		     fprintf(stdout, "\t%s\n", id);
		  else
		     fprintf(stdout, "\n");
               }
            }
            else
            {
               if (0 != strcmp(check_smiles, smiles))
               {
                  fprintf(stderr, "'%s' <> '%s'\n", check_smiles, smiles);
                  fprintf(stderr, "original line = %s\n", line);
                  // exit(1);
               }
            }
         }
// if (do_test) fprintf(stderr, "step 5\n");
         first_run=FALSE;
         MyFree(smiles);
         if (usmiles) MyFree(usmiles);
	 if (coordp)
	 {
	    MyFree((char *)coordp);
	    coordp = NULL;
	 }
// if (do_test) fprintf(stderr, "step 6\n");
         if (nruns == 2)   /* check for oscillations */
         {
            strcpy(buffer, check_smiles);
         }
         if (do_test && nruns < 3) goto repeat_do_test;
      }
// if (do_test) fprintf(stderr, "step 7\n");
      return (EXIT_SUCCESS);
   }

   else
   {
      explain(prognam, argc);
   }

   return (EXIT_SUCCESS);
}
