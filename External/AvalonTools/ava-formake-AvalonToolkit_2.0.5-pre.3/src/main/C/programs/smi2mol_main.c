/************************************************************************/
/*                                                                      */
/*    File:           smi2mol_main.c                                    */
/*                                                                      */
/*    Purpose:        This file implements the main executables that    */
/*                    convert SMILES and related formats to MOL file    */
/*                    and back.                                         */
/*                                                                      */
/*    History:        14-Feb-2020     Extracted from smi2mol.c          */
/*                                                                      */
/************************************************************************/
#include <stdlib.h>

#include <string.h>
#include <ctype.h>

#include "reaccs.h"
#include "utilities.h"
// #include "graph.h"
// #include "layout.h"
// #include "denormal.h"
// #include "stereo.h"
// #include "symbol_lists.h"
// #include "rtutils.h"
#include "perceive.h"
#include "smi2mol.h"
#include "shortcut.h"
#include "smi2mol.h"

#include "reaccsio.h"
#include "didepict.h"

#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

#define MXLN        2000
#define MXBUF        3000

static void explain(char * prognam, int argc);

FILE *log_file;

int main(int argc, char *argv[])
{
   FILE *fp;
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   struct reaccs_reaction_t *rp;
   char *smiles;
   char key[255];
   char *coordp;
   char line[MXLN];
   char buffer[MXBUF];
   char *cp, *cp1;
   char prognam[256];
   int nsmi, i, n, imol;
   char *data_columns[100];

   int show_help       = FALSE;
   int show_molname    = FALSE;
   int show_header     = FALSE;
   int show_coordinates= TRUE;
   int remove_mapping  = FALSE;
#define FP_SIZE 64
   char fp_buffer[FP_SIZE];
   int add_hex_fp  = FALSE;

   int more_flags = 0;
   int shortcut_flags = 0;
   int unknown_flag = FALSE;

   log_file = stderr;

   /* Read options */
   while (argc > 1  &&  argv[1][0] == '-' && argv[1][1] != '\0' )
   {
      if (0 == strcmp(argv[1], "-h"))
         show_help    = TRUE;
      else if (0 == strcmp(argv[1], "-x"))
         add_hex_fp = TRUE;
      else if (0 == strcmp(argv[1], "-r"))
         remove_mapping = TRUE;
      else if (0 == strcmp(argv[1], "-m"))
         show_molname = TRUE;
      else if (0 == strcmp(argv[1], "-s"))
         show_header  = TRUE;
      else if (0 == strcmp(argv[1], "-n"))
         show_coordinates  = FALSE;
      else if (0 == strcmp(argv[1], "-z"))
         more_flags   |= USE_Z;
      else if (0 == strcmp(argv[1], "-d"))
         more_flags   |= DY_AROMATICITY;
      else if (0 == strcmp(argv[1], "-shortcut"))
      {
         more_flags   |= PERCEIVE_SHORTCUTS;
         shortcut_flags   |= (AMINO_ACIDS|STANDARD_SHORTCUT);
      }
      else if (0 == strcmp(argv[1], "-sc_all"))
      {
         shortcut_flags   |= (EXTENDED_SHORTCUT | NON_STANDARD_SHORTCUT | N_METHYL_SHORTCUT | CATCH_ALL);
      }
      else if (0 == strcmp(argv[1], "-sc_extended"))
         shortcut_flags   |= EXTENDED_SHORTCUT;
      else if (0 == strcmp(argv[1], "-sc_non_standard"))
         shortcut_flags   |= NON_STANDARD_SHORTCUT;
      else if (0 == strcmp(argv[1], "-sc_n_methyl"))
         shortcut_flags   |= N_METHYL_SHORTCUT;
      else if (0 == strcmp(argv[1], "-sc_catch_all"))
         shortcut_flags   |= CATCH_ALL;
      else
         unknown_flag = TRUE;
      for (i=2; i<argc; i++) argv[i-1] = argv[i];
      argc--;
   }

   coordp = NULL;
   strncpy(prognam, argv[0], 255);
   for (cp=prognam; (*cp); cp++)
      (*cp) = tolower(*cp);


   if (show_help ||  unknown_flag)
   {
      explain(prognam, argc);
      return (EXIT_FAILURE);
   }


   /* convert first argument to MOL file */
   if (argc == 3 && (strstr(prognam, "smi2mol") || strstr(prognam, "sma2mol")))
   {
      if (strstr(prognam, "sma2mol"))
         mp = SMIToMOL(argv[1], DO_LAYOUT | EXPECT_SMARTS);
      else
         mp = SMIToMOL(argv[1], DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
      if (!mp) return (EXIT_FAILURE);

      fp = RedirFopen(argv[2],"w");
      /* patch for MOL file output to ISISDraw */
      for (i=0; i<mp->n_atoms; i++)
      {
         mp->atom_array[i].x *= 0.5;
         mp->atom_array[i].y *= 0.5;
      }
      PrintREACCSMolecule(fp, mp, "");
      FreeMolecule(mp);
// FORTIFY Fortify_ListAllMemory();
      return (EXIT_SUCCESS);
   }
   else if (argc == 3  &&  strstr(prognam, "smi2plt"))
   {
      fp = RedirFopen(argv[2],"w");
      if (!DepictSmilesToPlotFile(fp, argv[1], 0, 0)) return (EXIT_FAILURE);
// FORTIFY Fortify_ListAllMemory();
      return (EXIT_SUCCESS);
   }
   else if (argc == 2 && /* read stdin for reaction SMILES and write RDF to argv[1] */
            strstr(prognam, "smi2rdf"))
   {
      fp = RedirFopen(argv[1],"w");
      nsmi = 0;
      fprintf(fp, "$RDFILE 1\n");
      SetStripZeroes(TRUE);
      while (!feof(stdin))
      {
         nsmi++;
         fgets(line, MXLN-1, stdin);
         key[0] = buffer[0] = '\0';
         sscanf(line, "%s %s", buffer, key);
         if (feof(stdin)) break;

         rp = SMIToRXN(buffer);
         if (!rp) continue;

         strcpy(rp->program_name, "smi2rdf");
         if (key[0])
         {
            fprintf(fp, "$RFMT $RIREG %d\n", nsmi);
            strncpy(rp->name, key, MAXNAME); rp->name[MAXNAME] = '\0';
            PrintREACCSReaction(fp,rp);
            fprintf(fp, "$DTYPE ROOT:ID\n");
            fprintf(fp, "$DATUM %s\n", key);
         }
         else
         {
            fprintf(fp, "$RFMT $RIREG %d\n", nsmi);
            PrintREACCSReaction(fp,rp);
            fprintf(fp, "$DTYPE RXN:KEY\n");
            fprintf(fp, "$DATUM %d\n", nsmi);
         }
         FreeReaction(rp);
      }
      RedirFclose(fp);
// FORTIFY Fortify_ListAllMemory();
      return (EXIT_SUCCESS);
   }
   else if (argc == 1 && strstr(prognam, "smi2smi"))    /* read SMILES from stdin and write transformed SMILES to stdout */
   {
      nsmi = 0;
      while (!feof(stdin))
      {
         fgets(line, MXLN-1, stdin);
         key[0] = buffer[0] = '\0';
         /* SMILES can be followed by an idetifiere (key) to be added as a data field */
         sscanf(line, "%s %s", buffer, key);
         if (feof(stdin)) break;
         nsmi++;
         mp = SMIToMOL(buffer, DROP_TRIVIAL_HYDROGENS|DO_LAYOUT);
         if (more_flags & PERCEIVE_SHORTCUTS)
         {
             InitShortcuts();
// PrintREACCSMolecule(stderr,mp,"BEFORE");
             mp = ApplyShortcuts(mp, shortcut_flags);
// PrintREACCSMolecule(stderr,mp,"AFTER");
         }
         cp = MOLToSMIExt(mp, ISOMERIC_SMILES, (int *)NULL, &coordp);
         if (coordp)
         {
            MyFree((char *)coordp);
            coordp = NULL;
         }
         if (mp != NULL)
         {
             FreeMolecule(mp);
         }
         if (cp == NULL) continue;
         fprintf(stdout, "%s\n", cp);
         MyFree(cp);
      }
   }
   else if (argc == 2 && /* read stdin for SMILES and write SDF to argv[1] */
            strstr(prognam, "smi2mol"))
   {
      fp = RedirFopen(argv[1],"w");
      nsmi = 0;
      SetStripZeroes(TRUE);
      while (!feof(stdin))
      {
         nsmi++;
         fgets(line, MXLN-1, stdin);
         key[0] = buffer[0] = '\0';
         /* SMILES can be followed by an idetifiere (key) to be added as a data field */
         sscanf(line, "%s %s", buffer, key);
         if (feof(stdin)) break;

         mp = SMIToMOL(buffer, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
         if (!mp) continue;

         if (key[0])
         {
            strncpy(mp->name, key, 80);
            strncpy(mp->comment, key, 80);
            PrintREACCSMolecule(fp,mp,"");
            fprintf(fp, "> <key>\n%s\n\n$$$$\n", key);
         }
         else
         {
            PrintREACCSMolecule(fp,mp,"");
            fprintf(fp, "$$$$\n");
         }
         FreeMolecule(mp);
      }
      RedirFclose(fp);
// FORTIFY Fortify_ListAllMemory();
      return (EXIT_SUCCESS);
   }
   else if (argc == 4  && strstr(prognam, "mol2smi"))
   {      /* convert MOL-File or SD-File to SMILES with a key */
      ffp = FortranOpen(argv[1],"r");
      fp = RedirFopen(argv[2],"w");
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      imol = 0;
      while (FORTRAN_NORMAL == ReadREACCSMolecule(ffp,mp,""))
      {
// for (i=0; i<mp->n_atoms; i++)
// if (mp->atom_array[i].charge != 0) fprintf(stderr, "atom %d has charge %d\n", i+1, mp->atom_array[i].charge);
         smiles = MOLToSMIExt(mp,
                              ISOMERIC_SMILES|more_flags,
                              (int *)NULL,
                              &coordp);
         while (ffp->status == FORTRAN_NORMAL  &&   /* Search '$$$$' line */
                0 != strcmp(ffp->buffer, "$$$$"))   /* and get key        */
         {
            // grab a data field to be used as a key
            if (ffp->buffer[0] == '>'  &&
                strstr(ToLower(ffp->buffer), ToLower(argv[3])))
            {
               GetBuffer(ffp);
               strcpy(key, ffp->buffer);
            }
            GetBuffer(ffp);
         }
         if (strcmp(argv[3],"-") == 0)
         {
            if (mp->name[0] > ' ') strcpy(key,mp->name);
            else                   sprintf(key, "%d", ++imol);
         }
         fprintf(fp,"%s\t%s\n",key,smiles);
         if (ffp->status != FORTRAN_NORMAL) break;
         GetBuffer(ffp);                            /* Skip '$$$$' line */
         FreeMolecule(mp);
         mp = TypeAlloc(1,struct reaccs_molecule_t);
         MyFree(smiles);
         if (coordp)
         {
            MyFree((char *)coordp);
            coordp = NULL;
         }
      }
      FreeMolecule(mp);
      RedirFclose(fp);
      if (ffp->status != FORTRAN_EOF)
         return (EXIT_FAILURE);
      else
      {
// FORTIFY Fortify_ListAllMemory();
         return (EXIT_SUCCESS);
      }
   }
   else if (strstr(prognam, "mol2tbl"))
   {  /* convert SD-File to SMILES with data as additional columns */
      if (show_header)
      {
         if (show_coordinates)
             fprintf(stdout, "SMILES\t2D");
         else
             fprintf(stdout, "SMILES");
         if (show_molname) fprintf(stdout, "\tMOLNAME");
         if (add_hex_fp) fprintf(stdout, "\tFP");
         for (i=1; i<argc; i++) fprintf(stdout, "\t%s", argv[i]);
         fprintf(stdout, "\n");
      }
      for (i=0; i<argc; i++) data_columns[i] = TypeAlloc(4001, char);

      ffp = FortranOpen("-","r");
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      while (FORTRAN_NORMAL == ReadREACCSMolecule(ffp,mp,""))
      {
          if (remove_mapping)
          {
              for (i=0; i<mp->n_atoms; i++)
                  mp->atom_array[i].mapping = NONE;
          }
         /* calculate SMILES and coordinates */
         smiles = MOLToSMIExt(mp,
                              ISOMERIC_SMILES|more_flags,
                              (int *)NULL,
                              &coordp);
         if (smiles) fprintf(stdout, "%s", smiles);
         if (show_coordinates)
         {
             if (coordp) fprintf(stdout, "\t%s", coordp);
             else        fprintf(stdout, "\t");
         }
         if (show_molname) fprintf(stdout, "\t%s", mp->name);
         if (add_hex_fp)
         {
            fprintf(stdout, "\t");
            SetFingerprintBits(mp, fp_buffer, FP_SIZE, USE_ALL_FEATURES|USE_NON_SSS_BITS, FALSE, USE_DY_AROMATICITY);
            SetFingerprintBits(mp, fp_buffer, FP_SIZE, USE_ALL_FEATURES|USE_NON_SSS_BITS, FALSE, ACCUMULATE_BITS);
            for (i=0; i<FP_SIZE; i++)
               fprintf(stdout,"%02X",0xFF&fp_buffer[i]);
         }

         /* grab data columns */
         for (i=1; i<argc; i++)
            if (argv[i][0] == '{')        /* common value */
            {
               cp = strstr(argv[i], "}");
               if (cp != NULL)
               {
                  n = (cp-argv[i])-1;
                  strncpy(data_columns[i], argv[i]+1, n);
                  data_columns[i][n] = '\0';
               }
               else
                  data_columns[i][0] = '\0';
            }
            else
               data_columns[i][0] = '\0';
         while (ffp->status == FORTRAN_NORMAL  &&   /* Search '$$$$' line */
                0 != strcmp(ffp->buffer, "$$$$"))   /* and get key        */
         {
            if (ffp->buffer[0] == '>'  &&
                (cp = strstr(ffp->buffer, "<")) != NULL  &&
                (cp1 = strstr(cp+1, ">")) != NULL) /* This is a data section */
            {
               strncpy(key,cp+1,(cp1-cp)-1); key[(cp1-cp)-1] = '\0';
               for (i=1; i<argc; i++)
                  if (0 == strcmp(ToLower(key), ToLower(argv[i])))
                  {
                     GetBuffer(ffp);
                     if (ffp->status == FORTRAN_NORMAL)
                        strncpy(data_columns[i], ffp->buffer, 4000);
                     break;
                  }
            }
            GetBuffer(ffp);
         }
         for (i=1; i<argc; i++)
            fprintf(stdout,"\t%s",data_columns[i]);
         fprintf(stdout, "\n");
         if (ffp->status != FORTRAN_NORMAL) break;
         GetBuffer(ffp);                            /* Skip '$$$$' line */
         FreeMolecule(mp);
         mp = TypeAlloc(1,struct reaccs_molecule_t);
         MyFree(smiles);
         if (coordp)
         {
            MyFree((char *)coordp);
            coordp = NULL;
         }
      }
      FreeMolecule(mp);
      if (ffp->status != FORTRAN_EOF)
         return (EXIT_FAILURE);
      else
      {
// FORTIFY Fortify_ListAllMemory();
         return (EXIT_SUCCESS);
      }
   }
   else if (argc == 3 && strstr(prognam, "rdf2smi"))
   {      /* convert RDF-File to SMILES w/o key */
      ffp = FortranOpen(argv[1],"r");
      fp = RedirFopen(argv[2],"w");
      fprintf(stderr,"reading %s\n",argv[1]);
      fprintf(stderr,"writing %s\n",argv[2]);
      while ((char *)(rp = ReadREACCSReaction(ffp)) != NULL)
      {
         for (mp=rp->reactants; mp != NULL; mp=mp->next)
         {
            for (i=0; i<mp->n_atoms; i++)
                mp->atom_array[i].mapping = NONE;
            smiles = MOLToSMIExt(mp,
                                 ISOMERIC_SMILES|more_flags,
                                 (int *)NULL,
                                 &coordp);
            fprintf(fp,"%s",smiles);
            fprintf(fp,"\t");
            if (coordp)
            {
               MyFree((char *)coordp);
               coordp = NULL;
            }
            if (smiles != NULL) MyFree((char *)smiles);
         }
         fprintf(fp,">>");
         for (mp=rp->products; mp != NULL; mp=mp->next)
         {
            for (i=0; i<mp->n_atoms; i++)
                mp->atom_array[i].mapping = NONE;
            smiles = MOLToSMIExt(mp,
                                 ISOMERIC_SMILES|more_flags,
                                 (int *)NULL,
                                 &coordp);
            fprintf(fp,"\t");
            fprintf(fp,"%s",smiles);
            if (coordp)
            {
               MyFree((char *)coordp);
               coordp = NULL;
            }
            if (smiles != NULL) MyFree((char *)smiles);
         }
         fprintf(fp,"\n");
         FreeReaction(rp);
      }
      RedirFclose(fp);
      fprintf(stderr, "status = %d\n", ffp->status);
      if (ffp->status != FORTRAN_EOF)
         return (EXIT_FAILURE);
      else
      {
// FORTIFY Fortify_ListAllMemory();
         return (EXIT_SUCCESS);
      }
   }
   else if (argc == 3 && strstr(prognam, "mol2smi"))
   {      /* convert MOL-File or SD-File to SMILES w/o key */
      ffp = FortranOpen(argv[1],"r");
      fp = RedirFopen(argv[2],"w");
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      while (FORTRAN_NORMAL == ReadREACCSMolecule(ffp,mp,""))
      {
         smiles = MOLToSMIExt(mp,
                              ISOMERIC_SMILES|more_flags,
                              (int *)NULL,
                              &coordp);
         fprintf(fp,"%s\n",smiles);
         while (ffp->status == FORTRAN_NORMAL  &&   /* Search '$$$$' line */
                0 != strcmp(ffp->buffer, "$$$$"))
            GetBuffer(ffp);
         if (ffp->status != FORTRAN_NORMAL) break;
         GetBuffer(ffp);                            /* Skip '$$$$' line */
         if (coordp)
         {
            MyFree((char *)coordp);
            coordp = NULL;
         }
      }
      RedirFclose(fp);
      if (ffp->status != FORTRAN_EOF)
         return (EXIT_FAILURE);
      else
      {
// FORTIFY Fortify_ListAllMemory();
         return (EXIT_SUCCESS);
      }
   }
   else if (argc == 2 &&
            (strstr(prognam,"mol2smi") || strstr(prognam, "mol2sma")))
   {      /* convert MOL-File 0r SD-File to SMILES and write to stdout */
       if (more_flags & PERCEIVE_SHORTCUTS) InitShortcuts();

       ffp = FortranOpen(argv[1],"r");
       mp = TypeAlloc(1,struct reaccs_molecule_t);
       while (FORTRAN_NORMAL == ReadREACCSMolecule(ffp,mp,""))
       {
          if (strstr(prognam, "mol2sma"))
              smiles = MOLToSMI(mp, ISOMERIC_SMILES | SMARTS_PERCEPTION);
          else
          {
              if (more_flags & PERCEIVE_SHORTCUTS)
                 mp = ApplyShortcuts(mp, shortcut_flags);
// PrintREACCSMolecule(stderr, mp, "");
              smiles = MOLToSMIExt(mp, ISOMERIC_SMILES|more_flags, (int *)NULL, &coordp);
          }
          printf("%s\n",smiles);
          while (ffp->status == FORTRAN_NORMAL  &&   /* Search '$$$$' line */
                 0 != strcmp(ffp->buffer, "$$$$"))
             GetBuffer(ffp);
          if (coordp)
          {
             MyFree((char *)coordp);
             coordp = NULL;
          }
          if (smiles) MyFree((char *)smiles);
          if (ffp->status != FORTRAN_NORMAL) break;
          GetBuffer(ffp);                            /* Skip '$$$$' line */
       }
       if (ffp->status != FORTRAN_EOF)
          return (EXIT_FAILURE);
       else
       {
          if (mp->n_atoms > 0) FreeMolecule(mp);
          else                 MyFree((char *)mp);
// FORTIFY Fortify_ListAllMemory();
          return (EXIT_SUCCESS);
       }
   }
   else
   {
      explain(prognam, argc);
   }

// FORTIFY Fortify_ListAllMemory();
   return (EXIT_SUCCESS);
}

static void explain(char * prognam, int argc)
{
   if (strstr(prognam, "mol2smi"))
   {
      fprintf(stderr,"mol2smi used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: mol2smi [-d] [-shortcut] <SDF-file>\n");
      fprintf(stderr,"       or\n");
      fprintf(stderr,"       mol2smi [-d] <SDF-file> <smiles-table>\n");
      fprintf(stderr,"       or\n");
      fprintf(stderr,"       mol2smi [-d] <SDF-file> <smiles-table> <key name>\n");
      fprintf(stderr,"-d        performs Daylight aromaticity perception.\n");
      fprintf(stderr,"-shortcut perceives standard shortcut fragments and converts them to atext nodes.\n");
      exit (EXIT_FAILURE);
   }
   else if (strstr(prognam, "rdf2smi"))
   {
      fprintf(stderr,"rdf2smi used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: rdf2smi [-d] <RDF-file> <smiles-table>\n");
      exit (EXIT_FAILURE);
   }
   else if (strstr(prognam, "mol2tbl"))
   {
      fprintf(stderr,"mol2tbl used with %d arguments\n",argc-1);
      fprintf(stderr,
      "usage: mol2tbl [-h,-r,-m,-z,-d,-s,-x] key1 key2 ... <SDF-file >table-file\n");
      fprintf(stderr,"-s adds header row.\n");
      fprintf(stderr,"-r removes mapping.\n");
      fprintf(stderr,"-x add hex-coded fingerprint column.\n");
      fprintf(stderr,"-m shows MOLNAME.\n");
      fprintf(stderr,"-z uses z-coordinate for stereo.\n");
      fprintf(stderr,"-n suppress coordinates.\n");
      fprintf(stderr,"-d performs Daylight aromaticity perception.\n");
      exit (EXIT_FAILURE);
   }
   else if (strstr(prognam, "mol2sma"))
   {
      fprintf(stderr,"mol2sma used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: mol2sma <MOL file>\n");
      exit (EXIT_FAILURE);
   }
   else if (strstr(prognam,"smi2mol"))
   {
      fprintf(stderr,"smi2mol used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: smi2mol <smiles string> <MOL-file>\n");
      fprintf(stderr,"       or\n");
      fprintf(stderr,"       smi2mol <output SDF-file>\n");
      exit (EXIT_FAILURE);
   }
   else if (strstr(prognam,"smi2rdf"))
   {
      fprintf(stderr,"smi2rdf used with %d arguments\n",argc-1);
      fprintf(stderr,"usage: smi2rdf <output RDF-file>\n");
      exit (EXIT_FAILURE);
   }
   else
   {  fprintf(stderr,"unknown program name %s\n",prognam);
      exit (EXIT_FAILURE);
   }
}
