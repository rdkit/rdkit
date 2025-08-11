/************************************************************************/
/*                                                                      */
/*  File:     struchk.c                                                 */
/*                                                                      */
/*  Purpose:  Driver program for checking structures from an SD-Files   */
/*            using a table of allowed augmented atoms and              */
/*            transformations.                                          */
/*                                                                      */
/*            The result of the checks is appended to each SD-File      */
/*            entry as a set of data fields.                            */
/*                                                                      */
//  Revisions: 1.0.1   25-Jun-92       Added version output 'usage'.    */
//                                                                      */
//            1.0.9   17-Jul-92       Started adding stereochemistry    */
//                                    checks to the code. In addition,  */
//                                    changing command line processing  */
//                                    to getopt().                      */
//                                                                      */
//            1.1.0   27-Jul-92       Completed enhancement requests    */
//                                    concerning new options, stereo,   */
//                                    and geometry errors.              */
//                                                                      */
//            1.9.9   11-Sep-92       Started reorganization of tool    */
//                                    set into one program with         */
//                                    several options or command file   */
//                                    input.                            */
//                                                                      */
//                    24-Sep-92       Continuing "Grand Unification"    */
//                                                                      */
//                    01-Oct-92       Added ForceStereo option          */
//                                                                      */
//            2.0.0   22-Feb-93       Changed user interface to         */
//                                    multiple letter options.          */
//								                                        */
//	     2.2.1   23-Apr-98        Prepared for Oracle call-outs         */
//                                                                      */
//	     2.2.2   20-Jan-10        Fixed an error message                */
//                                                                      */
//	     2.2.3   07-Jan-14        Allow quoted options in profile       */
//                                    Ignore more Sgroup lines on read  */
//                                                                      */
//	     2.2.4   12-Feb-14        Implemented option for tautomer       */
//                                    standardization.                  */
//                                                                      */
//	     2.2.5   04-May-22        Fixed reading of isotopic labels      */
//                                                                      */
/************************************************************************/

static char struchk_version[] = "2.2.5";

#ifdef __TURBOC__
#include <process.h>
unsigned _stklen    = 0xFEEE;
unsigned _ovrbuffer = 0x2000;
#endif

#include "local.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "forio.h"
#include "reaccs.h"
#include "utilities.h"
#include "reaccsio.h"
#include "pattern.h"
#include "symboltable.h"
#include "stereo.h"
#include "aacheck.h"
#include "fixcharges.h"
#include "patclean.h"
#include "ssmatch.h"

#include "struchk.h"

#define MAX_OPT  100

// FILE *log_file;
// FILE *aa_log;

int ApplyTautomer(struct reaccs_molecule_t *mp, struct reaccs_molecule_t *from_tautomer, struct reaccs_molecule_t *to_tautomer)
{
   int map1, map2, tmap1, tmap2;
   int i, j;
   struct reaccs_atom_t *ap, *aph;
   struct reaccs_bond_t *bp, *bph;
   int *atom_colors;
   ssmatch_t *match;

   match = SSMatch(mp, from_tautomer, SINGLE_MATCH);
   if (IsNULL(match)) return (FALSE);
// fprintf(stderr,"found match for from_tautomer with %d atoms\n", from_tautomer->n_atoms);
   // We do have a match here so start serious work
   // Allocate atom color storage and remember atom colors
   atom_colors = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++) atom_colors[i] = mp->atom_array[i].color;
   // set mapping numbers of matching source pattern as atom color
   for (i=0; i<match->n_match; i++)
      mp->atom_array[match->match_atoms[i]].color = from_tautomer->atom_array[i].mapping;
   // scan for completely mapped bonds and replace bond order with mapped bond from to_tautomer
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      map1 = mp->atom_array[bp->atoms[0]-1].mapping;
      map2 = mp->atom_array[bp->atoms[1]-1].mapping;
      if (map1 <= 0  ||  map2 <= 0) continue;   // not completely mapped => no change
      // search for to_tautomer bond with same pair of mappings
      for (j=0, bph=to_tautomer->bond_array; j<to_tautomer->n_bonds; j++, bph++)
      {
         tmap1 = to_tautomer->atom_array[bph->atoms[0]-1].mapping;
         tmap2 = to_tautomer->atom_array[bph->atoms[1]-1].mapping;
         if (tmap1 <= 0  ||  tmap2 <= 0) continue;   // not completely mapped => no match
         if ((tmap1 == map1  &&  tmap2 == map2)  ||
             (tmap2 == map1  &&  tmap1 == map2))
         {
            if (bp->bond_type == bph->bond_type) break;  // no change
            bp->bond_type = bph->bond_type;
         }
      }
   }
   // apply charge/radical fixes if any
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->mapping == 0) continue;
      // search for to_tautomer atom with same mapping
      for (j=0, aph=to_tautomer->atom_array; j<to_tautomer->n_atoms; j++, aph++)
      {
         if (ap->mapping != aph->mapping) continue;
         if (ap->charge == aph->charge  &&  ap->radical == aph->radical) break;
         ap->charge = aph->charge; ap->radical = aph->radical;
      }
   }

   // Restore atom colors and free atom color storage
   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].color = atom_colors[i];
   MyFree((char *)atom_colors);
   FreeSSMatch(match);
   return (TRUE);
}

void ConvertGroupsToSGroups(struct reaccs_molecule_t *mp)
/*
 * Converts the line pairs defining ISIS groups to MACCS SGroup lines.
 * Currently, only lines which define fragment groups are translated.
 */
{
   struct prop_line_t head, *hp, *hhp, *tailp;
   int group_atom, linked_to;
   int nsgroups;

   nsgroups = 0;
   head.next = (struct prop_line_t *)NULL;
   tailp = &head;
   while (mp->prop_lines)
      if (STRING_BEGINS(mp->prop_lines->text,"G  ")  &&
          2 == sscanf(mp->prop_lines->text+strlen("G  "), "%d %d", &group_atom, &linked_to) &&
          linked_to == 0)
      {
         hp = mp->prop_lines->next;
         MyFree((char *)mp->prop_lines);
         mp->prop_lines = hp;
         mp->n_props--;
         if (hp)
         {
            mp->prop_lines = hp->next; mp->n_props--;
            nsgroups++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,"M  STY  1 %3d %-3s",nsgroups,"DAT");
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,"M  SLB  1 %3d %3d",nsgroups,nsgroups);
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,"M  SAL %3d %2d %3d",nsgroups,1,group_atom);
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,"M  SDT %3d %-3s",nsgroups,"COMPONENT");
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,
                    "M  SDD %3d     0.0000    0.0000    AR    ALL  1       7",
                    nsgroups);
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            hhp = TypeAlloc(1, struct prop_line_t);
            sprintf(hhp->text,"M  SED %3d %-3s",nsgroups,hp->text);
            tailp->next = hhp; tailp = hhp;
            mp->n_props++;

            MyFree((char *)hp);
         }
      }
      else     /* Copy untreated lines */
      {
         hp = mp->prop_lines; mp->prop_lines = hp->next;
         hp->next = (struct prop_line_t *)NULL; tailp->next = hp;
         tailp = hp;
      }

   mp->prop_lines = head.next;
}

int SqueezeSubstring(char *string, char *substring)
/*
 * Checks if string[] contains substring[]. If this
 * is the case, the substring is "squeezed" out and TRUE is
 * returned. If it was not the case, then the function returns FALSE.
 */
{
   char *cp, *cph;

   cp = strstr(string, substring);

   if (cp == NULL) return (FALSE);

   for (cph=cp+strlen(substring); *cph; cp++, cph++)
      (*cp) = (*cph);
   (*cp) = '\0';
   return (TRUE);
}

symbol_entry_t radical_symbols[] =
   {
      {TRIPLET, "^^^^"},
      {DOUBLET, "^^"},
      {SINGLET, ":"},
      {0, (char *)NULL}
   };

symbol_entry_t charge_symbols[] =
   {
      {15, "+15"}, {(-15), "-15"},
      {14, "+14"}, {(-14), "-14"},
      {13, "+13"}, {(-13), "-13"},
      {12, "+12"}, {(-12), "-12"},
      {11, "+11"}, {(-11), "-11"},
      {10, "+10"}, {(-10), "-10"},
      {9,  "+9"}, {(-9),  "-9"},
      {8,  "+8"}, {(-8),  "-8"},
      {7,  "+7"}, {(-7),  "-7"},
      {6,  "+6"}, {(-6),  "-6"},
      {5,  "+5"}, {(-5),  "-5"},
      {4,  "+4"}, {(-4),  "-4"},
      {3,  "+3"}, {(-3),  "-3"},
      {2,  "+2"}, {(-2),  "-2"},
      {1,  "+1"}, {(-1),  "-1"},
      {1,  "+"}, {(-1),  "-"},
      {0, (char *)NULL}
   };

int StringToAtomProperties(char *string,
                           char *symbol, int *charge, int *radical)
/*
 * Converts the ChemText atom text string[] to the corresponding
 * properties. The function returns TRUE if everything worked OK
 * and FALSE in case of error.
 */
{
   char buffer[256];
   symbol_entry_t *table;
   char *cp;

   strncpy(buffer,string,255); buffer[255] = '\0';

   for (table=radical_symbols; table->symbol_string; table++)
      if (SqueezeSubstring(buffer,table->symbol_string))
      {
	 (*radical) = table->symbol_id;
	 break;
      }

   for (table=charge_symbols; table->symbol_string; table++)
      if (SqueezeSubstring(buffer,table->symbol_string))
      {
         (*charge) = table->symbol_id;
	 break;
      }

   if (isalpha(buffer[0]) && isupper(buffer[0]))
   {
      for (cp = buffer+1; *cp; cp++)
	 if (!isalpha(*cp) || !islower(*cp)) return (FALSE);
      if (strlen(buffer) < 4)
      {
         strcpy(symbol, buffer);
         return (TRUE);
      }
      else
         return (FALSE);
   }
   else if (buffer[0] == 'r'  &&  '0' <= buffer[1]  &&  buffer[1] <= '9')
   {
      switch (buffer[1])
      {
         case '1':
            strcpy(symbol,"Ra");
            break;
         case '2':
            strcpy(symbol,"Rb");
            break;
         default:
            return (FALSE);
      }
      return (TRUE);
   }
   else
      return (FALSE);
}

int ConvertAtomAliases(struct reaccs_molecule_t *mp)
/*
 * Converts ChemText atom aliases to the corresponding
 * atom properties. It returns 1 if conversions went OK, 2 if no conversion were done, and 0 on failure.
 */
{
   struct prop_line_t head, *hp, *tailp;
   char symbol[4];
   int charge, radical;
   int atom;
   int result = 2;

   head.next = (struct prop_line_t *)NULL;
   tailp = &head;
   while (mp->prop_lines)    /* scan property lines for atom texts */
      if (STRING_BEGINS(mp->prop_lines->text,"A  ")  &&
          1 == sscanf(mp->prop_lines->text+strlen("A  "), "%d", &atom) &&
          atom > 0)
      {
         if (result > 1) result = 1;
         hp = mp->prop_lines->next;     /* remove atom text header */
         MyFree((char *)mp->prop_lines);
         mp->prop_lines = hp;
         mp->n_props--;

         if (hp)                        /* remove and interpret text */
         {
            mp->prop_lines = hp->next; mp->n_props--; hp->next = NULL;

            symbol[0] = '\0'; charge = 0; radical = NONE;
            if (!StringToAtomProperties(hp->text,
                                        symbol, &charge, &radical))
            {
               AddMsgToList("atom alias syntax error");
               AddMsgToList(hp->text);
               result = 0;
            }
            if (symbol[0] != '\0')
               strcpy(mp->atom_array[atom-1].atom_symbol,symbol);
            if (charge    != 0)
               mp->atom_array[atom-1].charge  = charge;
            if (radical   != NONE)
               mp->atom_array[atom-1].radical = radical;

            MyFree((char *)hp);
         }
         else
         {
            AddMsgToList("atom alias line missing");
            result = 0;
         }
      }
      else     /* Copy untreated lines */
      {
         hp = mp->prop_lines; mp->prop_lines = hp->next;
         hp->next = (struct prop_line_t *)NULL; tailp->next = hp;
         tailp = hp;
      }

   mp->prop_lines = head.next;
   return (result);
}

struct reaccs_molecule_t *StripSmallFragments(struct reaccs_molecule_t *mp, int *fragments_found)
/*
 * Splits the molecule *mp into its connected components and removes
 * all but the largest of them. *fragments_found is set TRUE if
 * fragments have been removed.
 */
{
   struct reaccs_molecule_t *result, *hp;

   *fragments_found = FALSE;
   mp->next = (struct reaccs_molecule_t *)NULL;
   mp = SplitMoleculeList(mp);
   result = mp; mp = mp->next; result->next = NULL;
   while (mp)
   {
      *fragments_found = TRUE;
      if (result->n_atoms < mp->n_atoms)
      {
         sprintf(msg_buffer,
                 "%10s    removed %d atoms",
                 result->name,
                 result->n_atoms);
         AddMsgToList(msg_buffer);
         FreeMolecule(result);
         result = mp; mp = mp->next; result->next = NULL;
      }
      else
      {
         hp = mp; mp = mp->next; hp->next = NULL;
         sprintf(msg_buffer,
                 "%10s    removed %d atoms",
                 result->name,
                 hp->n_atoms);
         AddMsgToList(msg_buffer);
         FreeMolecule(hp);
      }
   }

   return (result);
}

/*
 * Appends SDF data field records for molecular weight and formular
 * to data_list using the tags mw_tag and mf_tag, respectively.
 *
 * Note: The function returns the modified list. The actual use may require
 * that data_list be non-NULL.
 */
struct data_line_t *AddMWMF(struct data_line_t *data_list,
                            struct reaccs_molecule_t *mp,
                            char *mw_tag, char *mf_tag)
{
   struct data_line_t *dph;
   struct data_line_t *dphh;
   dphh = (struct data_line_t *)NULL;

   /* Add Molecular Weight */
   dph = TypeAlloc(1, struct data_line_t);
   strncpy(dph->data, "", MDL_MAXLINE);
   sprintf(dph->data, "%g\n", MolecularWeight(mp));
   dph->next = dphh; dphh = dph;
   dph = TypeAlloc(1, struct data_line_t);
   sprintf(dph->data, ">  <%s>", mw_tag);
   dph->next = dphh; dphh = dph;

   /* Add Molecular Formula */
   dph = TypeAlloc(1, struct data_line_t);
   MolecularFormula(mp, dph->data, MAXDATA+1-1);        // make sure the '\n' still fits into the buffer!
   strcat(dph->data, "\n");
   dph->next = dphh; dphh = dph;
   dph = TypeAlloc(1, struct data_line_t);
   sprintf(dph->data, ">  <%s>", mf_tag);
   dph->next = dphh; dphh = dph;

   if (data_list == NULL) return (dphh);

   /* find tail of list and append */
   for (dph = data_list; dph->next != NULL; dph=dph->next)
      ;
   dph->next = dphh;

   return (data_list);
}

void Cinderella(FILE *fp,
                struct reaccs_molecule_t *mp,
                struct data_line_t *data_list,
                int result_as_data,
                char *result)
/*
 * Writes the structure *mp to the file *fp, and appends the data lines
 * pointed * to by lp. If result_as_data is TRUE, it adds data lines with
 * the pending messages as the SD-entry <STRUCHK_MSG>, *result as
 * <STRUCHK_RES>, and the struchk version as <STRUCHK_VER>.
 */
{
   struct data_line_t *dph;

   PrintREACCSMolecule(fp,mp,"");
   for (dph = data_list; !IsNULL(dph); dph = dph->next)
      fprintf(fp,"%s\n",dph->data);
   if (result_as_data)
   {
      fprintf(fp, "> <STRUCHK_RES>\n%s\n\n",result);
      if (MsgsPending())
      {
         fprintf(fp,"> <STRUCHK_MSG>\n");
         PrintMsgs(fp);
         fprintf(fp,"\n");
      }
      fprintf(fp, "> <STRUCHK_VER>\n");
      fprintf(fp, "%s", struchk_version);
      fprintf(fp, " %s", aa_trans_version);
      fprintf(fp, " %s", aa_check_version);
      fprintf(fp, " %s", pKa_version);
      fprintf(fp, "\n\n");
   }
   fprintf(fp,"$$$$\n");
}

char *help_text[] =
{
   "The program is called with the following syntax:",
   " ",
   "     struchk {-option [parameter]}",
   "",
   "I. e., the command is followed by one letter option characters",
   "each optionally followed by a parameter. The options define",
   "the kind of processing and output desired.",
   "",
   "The following options are currently defined:",
   "",
   " -i <input file> Use this file for input. If no input file",
   "                 is given, struchk will print this message",
   "",
   " -o <output file>",
   "                 Use this file to output SD-file.",
   "                 Default output is stdout.",
   "",
   " -ov             Write version information to stdout",
   "",
   " -or             Show check result as data field.",
   "",
   " -b <bad file>   Write the \"bad\" structures to <bad file>.",
   "",
   " -br             Remove \"bad\" entries from output.",
   "",
   " -s <stereo output>",
   "                 Send every structure that has bee edited",
   "                 with respect to stereochemistry also to",
   "                 this file. Mainly used for EITHER bonds.",
   "",
   " -l <log file>   Write logging information to this file.",
   "",
   " -la <augmented atom log file>",
   "                 If augmented atoms are checked, any augmented",
   "                 atom found to be illegal is printed to this",
   "                 file. The format is like the input format for",
   "                 the check table file.",
   "",
   " -f <profile name>",
   "                 Use this file as option profile.",
   "",
   " -ft <tmp directory location>",
   "                 Use this directory to store temporary files.",
   "",
   " -ta <translation table file>",
   "                 Use augmented atom translations from this file.",
   "",
   " -tc <charge table file>",
   "                 Read pKa increments from this file. Don't",
   "                 recharge if this switch is not set. Structures",
   "                 which can not reliably be recharged are",
   "                 marked as BAD.",
   "",
   " -tl <acidity limit>",
   "                 Use this pKa limit for atom considered to",
   "                 be available for deprotonation.",
   "",
   " -tm             Split off minor fragments and keep only largest one",
   "",
   " -tq <desired charge>",
   "                 Sets the charge to be reached through",
   "                 deprotonation to <total charge>. This can be used",
   "                 to handle the registration of ionic species",
   "                 and non-stoichiometric representations.",
   "",
   " -tt <tautomer transformation file>",
   "                 RD-File containing tautomer transformations to be applied.",
   "",
   " -pc <pattern file>",
   "                 Use this file to read templates for",
   "                 pattern oriented cleaning.",
   "",
   " -pr <rotate pattern file>",
   "                 Use this file to read templates for",
   "                 pattern oriented rotation alignment.",
   "",
   " -ps <stereo pattern file>",
   "                 Force stereochemistry of matching templates.",
   "",
   " -ca <check table file>",
   "                 Check if all augmented atoms of the structures",
   "                 match one of the descriptions in the file.",
   "                 If one does not match, the structure is",
   "                 marked as BAD.",
   "",
   " -cc             Check for collisions of atoms with other",
   "                 atoms or bonds.",
   "",
   " -cl             Set distance limit for collision detection",
   "                 in percent of the average bond length.",
   "",
   " -cn <size>      Classify molecules with more than <size> atoms or",
   "                 bonds in the resulting structure as BAD.",
   "",
   " -cs             Check stereoconventions.",
   "",
   " -da             Convert atom text strings to properties.",
   "",
   " -dg             Translate ISIS groups to S-Groups.",
   "",
   " -dk <key>       Use the MACCS SD-file field <key> as the",
   "                 identification of the molecule during processing.",
   "",
   " -ds             Convert CPSS STEXT to data fields",
   "",
   " -dw             Squeeze whitespace out of idendtifiers",
   "",
   " -dz             Strip most of the trailing zeroes",
   "",
   " -sc             Remove CPSS properties ",
   "",
   " -sv             Enforce default valence if possible",
   "",
};

static int dontExit = FALSE;
static void MessageExit(char *format, char *program, char *parameter)
{
   char temp_name[1024];
   int i;

   fprintf(stderr, format, FileToCommandName(temp_name, program), parameter);
   if (!parameter)
   {
      for (i=0; i<sizeof(help_text)/sizeof(char *); i++)
         fprintf(stdout,"%s\n",help_text[i]);
   }
   if (!dontExit) exit (EXIT_FAILURE);
}

struct option_t
   {
      char *optstring;  /* option name string */
      int  has_par;     /* TRUE if option requires a parameter */
   }
   defined_options[] =
   {
      {"i",  TRUE},  /* input file */
      {"o",  TRUE},  /* output file */
      {"ov", FALSE}, /* output version */
      {"or", FALSE}, /* output results as data fields */
      {"b",  TRUE},  /* result file of bad structures */
      {"br", FALSE}, /* remove bad structure from output */
      {"s",  TRUE},  /* backup file for stereo edited structures */
      {"l",  TRUE},  /* log file */
      {"la", TRUE},  /* augmented atom check log file */
      {"f",  TRUE},  /* option source file */
      {"ft", TRUE},  /* name of temporary directory */
      {"ta", TRUE},  /* augmented atom translation table file */
      {"tc", TRUE},  /* deprotonation table file */
      {"tl", TRUE},  /* acidity limit */
      {"tm", FALSE}, /* remove all but largest fragment */
      {"tq", TRUE},  /* desired charge */
      {"tt", TRUE},  /* tautomerize patterns */
      {"pc", TRUE},  /* cleaning pattern file */
      {"pr", TRUE},  /* rotation pattern file */
      {"ps", TRUE},  /* pattern file for forced stereochemistry */
      {"ca", TRUE},  /* file of legal augmented atoms */
      {"cc", FALSE}, /* check for atom collisions */
      {"cl", TRUE},  /* distance limit for collision check in percent of average bond length */
      {"cn", TRUE},  /* limit to atom and bond count for good structures */
      {"cs", FALSE}, /* check MACCS stereoconventions */
      {"da", FALSE}, /* convert atom text strings to properties */
      {"dg", FALSE}, /* translate ISIS groups to S-groups */
      {"dk", TRUE},  /* data field to be used as MOLNAME */
      {"ds", FALSE}, /* Convert CPSS STEXT to data fields */
      {"dw", FALSE}, /* squeeze white space out of MOLNAME */
      {"dz", FALSE}, /* strip (most) trailing zeroes */
      {(char *)NULL, 0} /* end of table */
   };

int   opterr;
char *optstring;
char *optarg;
static int  optind = 1;

int GetOption(FILE *profile,
              int argc, char *argv[],
              struct option_t *options)
/*
 * Operates similar to the Unix getopt() function on the argc/argv
 * parameters, but it first reads options from profile. The profile
 * consists of lines, one line per option, with the following syntax:
 *
 *    <option string> [<parameter>] [<comment>]
 *
 * The function sets the global variable optarg just as getopt would.
 * It returns the first character of the option string (to be used
 * in a switch statement) and stores a pointer to the total option string
 * in optstring.
 *
 * When profile == NULL or end of file is reached, then the command
 * line options are processed.
 */
{
   static char option[20];
   static char parameter[1024];

   struct option_t *op;
   char buffer[1024];
   char *cp;
   int ipar;

   optstring = NULL;
   parameter[0] = '\0'; optarg = parameter;
   while (profile && !feof(profile) && fgets(buffer, 1023, profile))
   {
      if (0 >= sscanf(buffer,"%19s",option)) continue; /* skip blank lines */
      if (option[0] == '#') continue;                  /* skip comments */

      for (op = options; op->optstring; op++) /* option defined? */
         if (0 == strcmp(op->optstring, option))
            break;
      if (op->optstring == NULL)      /* undefined option */
      {
         optstring = "?";
         return ('?');
      }

      parameter[0] = '\0';
      if (op->has_par)
      {
         if (2 != sscanf(buffer,"%19s %1023s",option,parameter))
            optarg = NULL;
         else
         {
            if (parameter[0] != '"')
               optarg = parameter;
            else
            { // quoted parameters need special treatment
               // skip until character after "
               for (cp=buffer+strlen(option); (*cp) && cp[0] != '"'; cp++)
                  ;
               cp++;
               // copy until next " or EOL
               for (ipar=0; (*cp) && cp[0] != '"' && ipar < 1023; cp++)
                  parameter[ipar++] = (*cp);
               parameter[ipar] = '\0';
               optarg = parameter;
            }
         }
         optstring = option;
         return(option[0]);
      }
      else
      {
         optarg = (char *)NULL;
         optstring = option;
         return (option[0]);
      }
   }

   if (optind >= argc) return (-1);

   strcpy(option, argv[optind]+1);
   optstring = option; optind++;
   for (op = options; op->optstring; op++)  /* option defined? */
      if (0 == strcmp(op->optstring, option))
         break;
   if (op->optstring == (char *)NULL)            /* undefined option */
      return ('?');

   parameter[0] = '\0';
   if (op->has_par)
   {
      if (optind >= argc) return ('?');
      strcpy(parameter, argv[optind]);
      optarg = parameter; optind++;
   }
   else
      optarg = (char *)NULL;


   return (option[0]);
}

struct data_line_t *ConvertSTEXTToData(struct reaccs_molecule_t *mp,
                                       struct data_line_t *dp)
/*
 * Converts the STEXT lines in *mp to data lines and prepends them
 * to the existing list dp. It frees the space allocated to the
 * STEXT entries.
 */
{
   struct stext_line_t *sl;
   struct data_line_t *dph;

   if (mp->stext_lines)
   {
      dph = TypeAlloc(1, struct data_line_t);
      strncpy(dph->data, "", MDL_MAXLINE);
      dph->next = dp; dp = dph;
      while (mp->stext_lines)
      {
         dph = TypeAlloc(1, struct data_line_t);
         strncpy(dph->data, mp->stext_lines->text, MDL_MAXLINE);
         dph->next = dp; dp = dph;
         sl = mp->stext_lines;
         mp->stext_lines = sl->next;
         MyFree((char *)sl);
      }
      dph = TypeAlloc(1, struct data_line_t);
      strncpy(dph->data, ">  <CPSS_STEXT>", MDL_MAXLINE);
      dph->next = dp; dp = dph;
   }

   return (dp);
}


/***********************************************************************
   Own variables of this module. They will be initialized by InitCheckMol()
   and accessed by CheckMol.
************************************************************************/
/* -i <input file> */
static Fortran_FILE *fp = NULL;   /* input file */

/* -o <output file> */
static FILE *outfile = NULL;

/* -ov */
static int print_versions = FALSE;

/* -or */
static int result_as_data = FALSE;

/* -b <bad file> */
static FILE *bad_file = NULL;   /* copy of bad structures */

/* -br */
static int remove_bad_mols = FALSE;

/* -s <stereo output> */
static FILE *stereo_file = NULL;/* copy of str. with funny stereochemistry */

/* -f <profile> */
static FILE *profile = NULL;
static char *profile_name;

/* -ta <translation table file> */
static int ntrans = 0;
static aa_pair *trans_pairs = NULL;
static FILE *check_trans_file;

/* -tc <charge table file> */
static Fortran_FILE *charge_file = NULL; /* file with charge increments */
static int charges_read = FALSE;
static int charge_bad   = FALSE;
static int ndeprot, nrefine;
// Code for pKa optimization
// static char molform[255];
//

/* -tl <acidity limit> */
static double acidity_limit = 0.0;

/* -tm */
static int remove_minor_fragments = FALSE;

/* -tq <desired charge> */
static int desired_charge = 0;

/* -pc <pattern file> */
#define MAXPAT     10
static Fortran_FILE *pattern_file = NULL; /* file with clean patterns */
static struct reaccs_molecule_t *patterns[MAXPAT];
static int npat = 0;

/* -pr <rotate pattern file> */
#define MAXRPAT     10
static Fortran_FILE *rotate_pattern_file = NULL; /* file with rotate patterns */
static struct reaccs_molecule_t *rotate_patterns[MAXPAT];
static int nrpat = 0;

/* -ps */
#define MAXSTEREOPAT 10
static Fortran_FILE *stereo_pattern_file = NULL;
static struct reaccs_molecule_t *stereo_patterns[MAXSTEREOPAT];
static int nstereopat = 0;
static int stereo_bad = 0;

/* -tt <tautomere RD file> */
#define MAXTAUTOMERS  10
static Fortran_FILE *tautomer_file = NULL;   /* tautomer file */
static struct reaccs_molecule_t *from_tautomer[MAXTAUTOMERS];
static struct reaccs_molecule_t *to_tautomer[MAXTAUTOMERS];
static int ntautomers = 0;
static int ntautomerized = 0;

/* -ca <check file> */
static int ngood = 0;
static augmented_atom_t *good_atoms = NULL;

/* -cc */
static int check_collision = FALSE;

/* -cl <collision limit in %> */
static int percent = 15;

/* -cn <size> */
static int max_mol_size = 255;

/* -cs */
static int check_stereo = FALSE;
static int stereo_result;

/* -da */
static int convert_atom_texts = FALSE;

/* -dg */
static int groups_to_sgroups = FALSE;

/* -dk <key> */
static char key[80];           /* molecule idetification */

/* -ds */
static int convert_stext = FALSE;

/* -dw */
static int squeeze_identifiers = FALSE;

/* -dz */
static int strip_zeroes = FALSE;

/* -l */
FILE *log_file = (FILE *)NULL; /* global file */

/* -la */
FILE *aa_log = (FILE *)NULL; /* global file */

/* -ft */
extern char *tmp_dir_name; /* name for a tmp directory */

void ClearParameters()
{
   int i, j;

// fprintf(stderr,"ClearParameters(1)\n");
   if (!IsNULL(fp)) FortranClose(fp);
   fp = NULL;
   if (outfile != stdout)
   {
      if (!IsNULL(outfile)) fclose(outfile);
   }
   outfile = (FILE *)NULL;
// fprintf(stderr,"ClearParameters(2)\n");
   print_versions = FALSE;
   result_as_data = FALSE;
// fprintf(stderr,"ClearParameters(3)\n");
   if (!IsNULL(bad_file))
	   fclose(bad_file);
   bad_file = (FILE *)NULL;
   remove_bad_mols = FALSE;
   if (!IsNULL(stereo_file))
	   fclose(stereo_file);
   stereo_file = (FILE *)NULL;
   if (!IsNULL(profile))
	   fclose(profile);
   profile = (FILE *)NULL;
// fprintf(stderr,"ClearParameters(4)\n");
   profile_name = NULL;
   if (!IsNULL(trans_pairs))
   {
      for (i=0; i<ntrans; i++)
      {
         if (!IsNULL(trans_pairs[i]))
         {
            MyFree((char *)trans_pairs[i][0].atom_symbol);
            MyFree((char *)trans_pairs[i][0].short_name);
            for (j=0; j<trans_pairs[i][0].n_ligands; j++)
               MyFree((char *)trans_pairs[i][0].ligands[j].atom_symbol);
            MyFree((char *)trans_pairs[i][1].atom_symbol);
            MyFree((char *)trans_pairs[i][1].short_name);
            for (j=0; j<trans_pairs[i][1].n_ligands; j++)
               MyFree((char *)trans_pairs[i][1].ligands[j].atom_symbol);
         }
      }
// fprintf(stderr,"ClearParameters(4a)\n");
      MyFree((char *)trans_pairs); trans_pairs = NULL;
   }
// fprintf(stderr,"ClearParameters(5)\n");
   ntrans = 0;
   if (!IsNULL(check_trans_file)) fclose(check_trans_file);
   check_trans_file = (FILE *)NULL;
// fprintf(stderr,"ClearParameters(6)\n");
   if (!IsNULL(charge_file)) FortranClose(charge_file);
   charge_file = NULL;
// fprintf(stderr,"ClearParameters(7)\n");
   charge_file = NULL; /* file with charge increments */
   charges_read = FALSE;
   charge_bad   = FALSE;
   ndeprot = 0;
   nrefine = 0;
   acidity_limit = 0.0;
   remove_minor_fragments = FALSE;
   desired_charge = 0;
   if (!IsNULL(tautomer_file)) FortranClose(tautomer_file);
   tautomer_file = NULL;
   if (ntautomers != 0) {
      for (i=0; i<ntautomers; i++)
      {
         if (!IsNULL(from_tautomer[i])) FreeMolecule(from_tautomer[i]);
         from_tautomer[i] = (struct reaccs_molecule_t *)NULL;
         if (!IsNULL(to_tautomer[i])) FreeMolecule(to_tautomer[i]);
         to_tautomer[i] = (struct reaccs_molecule_t *)NULL;
      }
   }
   npat = 0;
   if (!IsNULL(pattern_file)) {
	   FortranClose(pattern_file);
   }
   pattern_file = NULL;
// fprintf(stderr,"ClearParameters(8)\n");
// static struct reaccs_molecule_t *patterns[MAXPAT];
// fprintf(stderr,"ClearParameters(9)\n");
   if (npat != 0)
      for (i=0; i<npat; i++)
      {
         if (!IsNULL(patterns[i])) FreeMolecule(patterns[i]);
         patterns[i] = (struct reaccs_molecule_t *)NULL;
      }
   npat = 0;
// fprintf(stderr,"ClearParameters(10)\n");
   if (!IsNULL(rotate_pattern_file)) FortranClose(rotate_pattern_file);
   rotate_pattern_file = NULL;
   rotate_pattern_file = NULL; /* file with rotate patterns */
   // static struct reaccs_molecule_t *rotate_patterns[MAXPAT];
// fprintf(stderr,"ClearParameters(11)\n");
   if (nrpat != 0)
      for (i=0; i<nrpat; i++)
      {
         if (!IsNULL(rotate_patterns[i])) FreeMolecule(rotate_patterns[i]);
         rotate_patterns[i] = (struct reaccs_molecule_t *)NULL;
      }
   nrpat = 0;
// fprintf(stderr,"ClearParameters(12)\n");
   if (!IsNULL(stereo_pattern_file)) FortranClose(stereo_pattern_file); stereo_pattern_file = NULL;
   stereo_pattern_file = NULL;
   // static struct reaccs_molecule_t *stereo_patterns[MAXSTEREOPAT];
   if (nstereopat != 0)
      for (i=0; i<nstereopat; i++)
      {
         if (!IsNULL(stereo_patterns[i])) FreeMolecule(stereo_patterns[i]);
         stereo_patterns[i] = (struct reaccs_molecule_t *)NULL;
      }
   nstereopat = 0;
   stereo_bad = 0;
   ntautomerized = 0;
// fprintf(stderr,"ClearParameters(13)\n");
   if (!IsNULL(good_atoms))
   {
      for (i=0; i<ngood; i++)
      {
         MyFree((char *)good_atoms[i].atom_symbol);
         MyFree((char *)good_atoms[i].short_name);
         for (j=0; j<good_atoms[i].n_ligands; j++)
            MyFree((char *)good_atoms[i].ligands[j].atom_symbol);
      }
      MyFree((char *)good_atoms);
   }
// fprintf(stderr,"ClearParameters(14)\n");
   ngood = 0;
   good_atoms = NULL;
   check_collision = FALSE;
   percent = 15;
   max_mol_size = 255;
   check_stereo = FALSE;
   stereo_result = 0;
   convert_atom_texts = FALSE;
   groups_to_sgroups = FALSE;
   key[0] = '\0';
   convert_stext = FALSE;
   squeeze_identifiers = FALSE;
   strip_zeroes = FALSE;
}

int Initialize(FILE *profile, int argc, char *argv[])
{
   int c;
   int i;
   struct reaccs_molecule_t *ssp;
   int errflg = 0;
   char *progname = "LIBRARY";

   if (argc <= 0) return 0;  // do not change anything when there is no input

   ClearParameters();

   if (argv[0][0] == '-') optind = 0;   // first string is already an option => start there
   else
   {
      optind = 1;   // first string is program name => skip it
      progname = argv[0];
   }
   key[0] = '\0';
   while ((-1) != (c = GetOption(profile,
                                 argc, argv,
                                defined_options)))
   { /* "i:t:c:q:p:s:mn:gk:wo:b:re:l:vdf:a:"))) */
// fprintf(stderr,"optstring = '%s', optarg = '%s'\n",optstring,optarg);
      switch (c)
      {
         case 'i':      /* input related switches */
            if (IsNULL(fp = FortranOpen(optarg,"r")))
	       MessageExit("%s: Could not open input file '%s'\n",
                           progname, optarg);
            break;
         case 'o':      /* output related switches */
	    if (strcmp(optstring,"o") == 0)     /* output file */
            {
               if (outfile != stdout)
               {
                  if (!IsNULL(outfile)) fclose(outfile);
               }
               outfile = (FILE *)NULL;
               if (IsNULL(outfile = file_open(optarg,"w")))
                  MessageExit("%s: Could not open output file '%s'\n",
                  progname, optarg);
            }
            else if (strcmp(optstring,"ov") == 0) /* output version */
               print_versions = TRUE;
            else if (strcmp(optstring,"or") == 0) /* results as data fields */
               result_as_data = TRUE;
            break;
         case 'b':      /* bad structure related switches */
            if (strcmp(optstring,"b") == 0)     /* output of "bad" str. */
            {
               if (!IsNULL(bad_file)) fclose(bad_file); bad_file = (FILE *)NULL;
               if (IsNULL(bad_file = file_open(optarg,"w")))
                  MessageExit("%s: Could not open file '%s'\n",
                              progname, optarg);
            }
            else if (strcmp(optstring,"br") == 0) /* remove bad mols */
               remove_bad_mols = TRUE;
            break;
         case 's':      /* stereo output file */
               if (!IsNULL(stereo_file)) fclose(stereo_file); stereo_file = (FILE *)NULL;
            if (IsNULL(stereo_file = file_open(optarg,"w")))
               MessageExit("%s: Could not open file '%s'\n", progname, optarg);
            break;
         case 'l':              /* log files */
            if (strcmp(optstring,"l") == 0)     /* general log file */
            {
               if (!IsNULL(log_file)) fclose(log_file); log_file = (FILE *)NULL;
               if (IsNULL(log_file = file_open(optarg,"w")))
                  MessageExit("%s: Could not open log file '%s'\n",
                              progname, optarg);
            }
            else if (strcmp(optstring,"la") == 0) /* augmented atom log */
            {
               if (!IsNULL(aa_log)) fclose(aa_log); aa_log = (FILE *)NULL;
               if (IsNULL(aa_log = file_open(optarg,"w")))
                  MessageExit("%s: Could not open log file '%s'\n",
                              progname, optarg);
            }
            break;
         case 'f':              /* option flag file */ /*AG*/
            if (strcmp(optstring,"f" ) == 0 )  /* option flag file */
            {
	       if (profile)
                  fclose(profile);
               if (IsNULL(profile = file_open(optarg,"r")))
                  MessageExit("%s: Could not open profile file '%s'\n",
                              progname, optarg);
            }
            else if (strcmp(optstring,"ft") == 0 ) /* tmp directory */
            {
	       if (*optarg == '\0' || optarg == NULL ) optarg = "/tmp";
               tmp_dir_name = TypeAlloc(strlen(optarg)+1,char);
               strcpy(tmp_dir_name,optarg);
               if (tmp_dir_name[strlen(tmp_dir_name)-1] == '/' )
                  tmp_dir_name[strlen(tmp_dir_name)-1] = '\0';
            }
	    break;
         case 't':              /* transformation options */
            if (strcmp(optstring,"ta") == 0)      /* augm. atom transf. */
            {
               if (IsNULL(check_trans_file = file_open(optarg,"r")))
                  MessageExit("%s: Could not open translation file '%s'\n",
                              progname, optarg);
               trans_pairs = ReadAAPairs(check_trans_file,&ntrans);
               if (IsNULL(trans_pairs))
                  MessageExit("%s: Error in translation file '%s'\n",
                              progname, optarg);
               fclose(check_trans_file); check_trans_file = (FILE *)NULL;
            }
            else if (strcmp(optstring,"tc") == 0) /* deprotonation */
            {
               if (IsNULL(charge_file = FortranOpen(optarg,"r")))
                  MessageExit("%s: Could not open charge file '%s'\n",
                              progname, optarg);
               if (!InitializeChargeDataTables(charge_file))
                  MessageExit("%s: Format error in charge file '%s'\n",
                              progname, optarg);
               charges_read = TRUE;
               FortranClose(charge_file); charge_file = NULL;

// Code for pKa optimization
// SetChargeLog(stdout);
// PrintChargeHeader();
//
            }
            else if (strcmp(optstring,"tq") == 0) /* desired charge */
               desired_charge = atoi(optarg);
            else if (strcmp(optstring,"tl") == 0) /* acidity limit */
            {
               acidity_limit = atof(optarg);
               if (acidity_limit > 0.0)
	       SetAcidityLimit(acidity_limit);
            }
            else if (strcmp(optstring,"tm") == 0) /* only largest fragment */
               remove_minor_fragments = TRUE;
            else if (strcmp(optstring,"tt") == 0) /* use tautomerization rules */
            {
               //
               ssp = TypeAlloc(1,struct reaccs_molecule_t);
               tautomer_file = FortranOpen(optarg,"r");
               if (IsNULL(tautomer_file = FortranOpen(optarg,"r")))
                  MessageExit("%s: Could not open tautomer rule file '%s'\n",
                              progname, optarg);
               ntautomers = 0;
               while (SearchString(tautomer_file, "$RXN", "##EOF##"))
               {
                  // skip to $MOL block
                  if (!SearchString(tautomer_file, "$MOL", "##EOF##")) break;
                  GetBuffer(tautomer_file);
                  // read pattern to be tautomerized
                  if (FORTRAN_NORMAL != ReadREACCSMolecule(tautomer_file,ssp,"")) break;
                  MakeHydrogensImplicit(ssp);
                  from_tautomer[ntautomers] = ssp;
// PrintREACCSMolecule(stderr,ssp,"from tautomer");
                  ssp = TypeAlloc(1, struct reaccs_molecule_t);
                  // skip to next $MOL block
                  if (!SearchString(tautomer_file, "$MOL", "##EOF##")) break;
                  GetBuffer(tautomer_file);
                  // read charge/bond order target pattern
                  if (FORTRAN_NORMAL != ReadREACCSMolecule(tautomer_file,ssp,"")) break;
                  MakeHydrogensImplicit(ssp);
// PrintREACCSMolecule(stderr,ssp,"to tautomer");
                  to_tautomer[ntautomers] = ssp;
                  ssp = TypeAlloc(1, struct reaccs_molecule_t);
                  ntautomers++;
break;
               }
               if (ntautomers == 0)
                  MessageExit("%s: Format error in tautomer pattern file '%s'\n",
			      progname, optarg);
               FreeMolecule(ssp);
               FortranClose(tautomer_file); tautomer_file = NULL;
            }
            break;
         case 'p':              /* pattern oriented options */
            if (strcmp(optstring,"pc") == 0)    /* pattern cleaning */
            {
               ssp = TypeAlloc(1,struct reaccs_molecule_t);
               pattern_file = FortranOpen(optarg,"r");
               if (IsNULL(pattern_file = FortranOpen(optarg,"r")))
                  MessageExit("%s: Could not open clean pattern file '%s'\n",
                              progname, optarg);
               npat = 0;
               while (FORTRAN_NORMAL == ReadREACCSMolecule(pattern_file,ssp,""))
               {
                  if (npat >= MAXPAT)
                  {
                     fprintf(stderr,"More than %d clean pattens\n",MAXPAT);
                     fprintf(stderr,"Using only the first %d of them\n",MAXPAT);
                     break;
                  }
                  while (pattern_file->status == FORTRAN_NORMAL)
                  {                                     /* Skip $$$$ line */
                     GetBuffer(pattern_file);
                     if (strncmp(pattern_file->buffer,"$$$$",4) == 0)
                     {
                        GetBuffer(pattern_file);
                        break;
                     }
                  }
                  MakeHydrogensImplicit(ssp);

// PrintREACCSMolecule(stderr,ssp,"");
// for (i=0; i<ssp->n_atoms; i++) ssp->atom_array[i].query_H_count = NONE;
                  patterns[npat] = ssp; npat++;
                  ssp = TypeAlloc(1, struct reaccs_molecule_t);
               }
               if (npat == 0)
                  MessageExit("%s: Format error in clean pattern file '%s'\n",
			      progname, optarg);
               FreeMolecule(ssp);
               FortranClose(pattern_file); pattern_file = NULL;
            }
            else if (strcmp(optstring,"pr") == 0)    /* pattern rotation */
            {
               ssp = TypeAlloc(1,struct reaccs_molecule_t);
               rotate_pattern_file = FortranOpen(optarg,"r");
               if (IsNULL(rotate_pattern_file = FortranOpen(optarg,"r")))
                  MessageExit("%s: Could not open rotate pattern file '%s'\n",
                              progname, optarg);
               nrpat = 0;
               while (FORTRAN_NORMAL ==
	              ReadREACCSMolecule(rotate_pattern_file,ssp,""))
               {
                  if (nrpat >= MAXRPAT)
                  {
                     fprintf(stderr,"More than %d rotate pattens\n",MAXRPAT);
                     fprintf(stderr,"Using only %d of them\n",MAXRPAT);
                     break;
                  }
                  while (rotate_pattern_file->status == FORTRAN_NORMAL)
                  {                                     /* Skip $$$$ line */
                     GetBuffer(rotate_pattern_file);
                     if (strncmp(rotate_pattern_file->buffer,"$$$$",4) == 0)
                     {
                        GetBuffer(rotate_pattern_file);
                        break;
                     }
                  }
                  MakeHydrogensImplicit(ssp);
                  for (i=0; i<ssp->n_atoms; i++)
                     ssp->atom_array[i].query_H_count = NONE;
                  rotate_patterns[nrpat] = ssp; nrpat++;
                  ssp = TypeAlloc(1, struct reaccs_molecule_t);
               }
               if (nrpat == 0)
                  MessageExit("%s: Format error in rotate pattern file '%s'\n",
			      progname, optarg);
               FreeMolecule(ssp);
               FortranClose(rotate_pattern_file); rotate_pattern_file = NULL;
            }
            else if (strcmp(optstring,"ps") == 0)       /* force stereo */
            {
               ssp = TypeAlloc(1,struct reaccs_molecule_t);
               stereo_pattern_file = FortranOpen(optarg,"r");
               if (IsNULL(stereo_pattern_file = FortranOpen(optarg,"r")))
                  MessageExit("%s: Could not open stereo pattern file '%s'\n",
                              progname, optarg);
               nstereopat = 0;
               while (FORTRAN_NORMAL ==
                      ReadREACCSMolecule(stereo_pattern_file,ssp,""))
               {
                  if (nstereopat >= MAXSTEREOPAT)
                  {
                     fprintf(stderr,
                             "More than %d stereo pattens\n",
                             MAXSTEREOPAT);
                     fprintf(stderr,
                             "Using only the first %d of them\n",
                             MAXSTEREOPAT);
                     break;
                  }
                  while (stereo_pattern_file->status == FORTRAN_NORMAL)
                  {                                     /* Skip $$$$ line */
                     GetBuffer(stereo_pattern_file);
                     if (strncmp(stereo_pattern_file->buffer,"$$$$",4) == 0)
                     {
                        GetBuffer(stereo_pattern_file);
                        break;
                     }
                  }
		  MakeHydrogensImplicit(ssp);
                  for (i=0; i<ssp->n_atoms; i++)
                     ssp->atom_array[i].query_H_count = NONE;
                  stereo_patterns[nstereopat] = ssp; nstereopat++;
                  ssp = TypeAlloc(1, struct reaccs_molecule_t);
               }
               if (nstereopat == 0)
                  MessageExit("%s: Format error in clean pattern file '%s'\n",
                              progname, optarg);
               FreeMolecule(ssp);
               FortranClose(stereo_pattern_file); stereo_pattern_file = NULL;
            }
            break;
         case 'c':              /* checking options */
            if (strcmp(optstring,"ca") == 0)      /* augm. atom checks */
            {
               if (IsNULL(check_trans_file = file_open(optarg,"r")))
                  MessageExit("%s: Could not open check file '%s'\n",
                              progname, optarg);
               good_atoms = ReadAugmentedAtoms(check_trans_file,&ngood);
               if (IsNULL(good_atoms))
                  MessageExit("%s: Error in augmented atom file '%s'\n",
                              progname, optarg);
               fclose(check_trans_file); check_trans_file = (FILE *)NULL;
            }
            else if (strcmp(optstring,"cc") == 0)  /* check for collisions */
               check_collision = TRUE;
            else if (strcmp(optstring,"cl") == 0) /* collision dist. limit */
               if (1 == sscanf(optarg,"%d",&percent))
                  SetCollisionLimit(percent);
               else
                  MessageExit("%s: argument of 'cl' (%s) not an integer\n",
                              progname, optarg);
            else if (strcmp(optstring,"cn") == 0) /* size limit */
               max_mol_size = atoi(optarg);
            else if (strcmp(optstring,"cs") == 0)  /* check stereoconv. */
               check_stereo = TRUE;
            break;
         case 'd':              /* data related options */
            if (strcmp(optstring,"da") == 0)      /* atom text to props. */
               convert_atom_texts = TRUE;
            else if (strcmp(optstring,"dg") == 0) /* ISIS- to S-groups */
            {
               groups_to_sgroups = TRUE;
            }
            else if (strcmp(optstring,"dk") == 0) /* MOLNAME field */
            {
               if (IsNULL(optarg))
                  errflg++;
               else
                  strcpy(key,optarg);
            }
            else if (strcmp(optstring,"ds") == 0) /* convert STEXTs */
            {
               convert_stext = TRUE;
            }
            else if (strcmp(optstring,"dw") == 0) /* squeeze white space */
            {
               squeeze_identifiers = TRUE;
            }
            else if (strcmp(optstring,"dz") == 0) /* strip zeroes */
            {
	       SetStripZeroes(TRUE);
               strip_zeroes = TRUE;
            }
            break;
         case '?':
            fprintf(stderr,"Illegal option '%s' with argument '%s'\n",
                    optstring, optarg);
            errflg++;
            break;
         default:
            fprintf(stderr,"Don't know option '%c'\n",c);
            errflg++;
            break;
      }
   }

   return (errflg);
}

/**
 * This function converts parameters from single, linefeed-separated
 * string representation to array of strings. This is an entry point
 * that is supposed to be called e.g. from the Oracle call-out glue code.
 *
 * Note: The opt string is modified and used to provide storage for '\0' terminated strings
 * collected as 'argument' strings.
 */
int InitCheckMol(char *opt)
{
   char *argv[MAX_OPT];   // up to MAX_OPT option strings
   char *cp;
   int argc=0;
   int flags;
   int hasCRLF, hasQuote, wordInLine, eolFound, eosFound;

   dontExit = TRUE;
   if (IsNULL(opt)) return 0;
   /* now, parse into argument possibly quoted strings */
   // check and remember if there might be comments after the optional second string on a line
   hasCRLF = strchr(opt, '\r') >= 0  ||  strchr(opt, '\n') >= 0;
// fprintf(stderr,"opt='%s', hasCRLF=%d\n",opt,hasCRLF);
   wordInLine = 0;
   do
   {
      hasQuote = FALSE; eolFound = FALSE;
      // skip until first non-trivial character
      while (isspace(opt[0]))
      {
          if (opt[0] == '\r' ||  opt[0] == '\n') eolFound = TRUE;
          opt++;
      }
      // skip and remember '"' if any
      if (opt[0] == '"')
      {
          hasQuote = TRUE;
          opt++;
      }
      if (hasQuote)
      {
          wordInLine++;
          for (cp = opt; cp[0] != '"'  &&  cp[0] != '\0'; cp++)
              ;
          eosFound =  (cp[0] == '\0');
          cp[0] = '\0';
       // fprintf(stderr, "1) argv[%d] = '%s'\n", argc, opt);
          argv[argc++] = opt;
          if (!eosFound) opt = cp+1;
      }
      else if (opt[0] != '\0')
      {
          wordInLine++;
          for (cp = opt; !isspace(cp[0])  &&  cp[0] != '\0'; cp++)
              ;
          eosFound =  (cp[0] == '\0');
          if (cp[0] == '\r' ||  cp[0] == '\n') eolFound = TRUE;

          cp[0] = '\0';
       // fprintf(stderr, "2) argv[%d] = '%s'\n", argc, opt);
          argv[argc++] = opt;
          if (!eosFound) opt = cp+1;
      }
      if (wordInLine == 2  &&  hasCRLF  &&  !eolFound)     // skip to after EOL
      {
          while (opt[0] != '\0'  &&  opt[0] != '\r' &&  opt[0] != '\n')
              opt++;
          if (opt[0] == '\r' ||  opt[0] == '\n')
              eolFound = TRUE;
          wordInLine = 0;
      }
      else if (eolFound)
      {
          wordInLine = 0;
      }
   } while (argc < MAX_OPT && opt[0]);
   flags = Initialize((FILE *)NULL, argc, argv);

   return (flags);
}

int _InitCheckMol_(char *opt)
{
   char **argv;   // up to MAX_OPT option strings
   char *arg, *tmp;
   int i, argc;
   int flags;
   int hasCRLF;

   if (IsNULL(opt)) return 0;
   /* now, parse into argument possibly quoted strings */
   // check and remember if there might be comments after the optional second string on a line
   hasCRLF = strchr(opt, '\r') >= 0  ||  strchr(opt, '\n') >= 0;
   argc = 0;
   /* count number of arguments */
   argc = 0;
   tmp = TypeAlloc(strlen(opt)+1, char);
   strcpy(tmp, opt);
   for (arg=strtok(tmp, " \t\r\n"); arg && (*arg); arg=strtok(NULL, " \t\r\n"))
      argc++;

   argv = TypeAlloc(argc+1, char *);
   strcpy(tmp, opt);
   /* create argument array */
   argv[0] = "";	/* GetOptions starts with 1 */
   for (i=1, arg=strtok(tmp, " \t\r\n"); arg && (*arg); arg=strtok(NULL, " \t\r\n"))
      argv[i++] = arg;

   flags = Initialize((FILE *)NULL, argc+1, argv);

   /* free memory */
   MyFree(tmp);
   // MyFree((char *)argv);
   return (flags);
}

int RunStruchk(struct reaccs_molecule_t **mpp, struct data_line_t *data_list)
{
   int i, j;
   int result = 0;
   int tmp;
   int fragments_found = 0;
   struct reaccs_molecule_t *ssp;
   struct reaccs_molecule_t *mp, *new_mp;
   struct data_line_t *new_data_list = (struct data_line_t *)NULL;
   struct data_line_t *dph;

   mp = (*mpp); // just used for easy writing

   if (mp->n_atoms > max_mol_size  ||  mp->n_bonds > max_mol_size  ||
       mp->n_atoms > MAXATOMS      ||  mp->n_bonds > MAXATOMS)
   {
      snprintf(msg_buffer, MAXMSG,
	      "%10s    : more than %d atoms or bonds",
	      mp->name, MAXATOMS>max_mol_size?max_mol_size:MAXATOMS);
      AddMsgToList(msg_buffer);
      result |= SIZE_CHECK_FAILED;
   }

   if ((result & SIZE_CHECK_FAILED) == 0)
   {
      for (i = 0; i < mp->n_bonds; ++i) {
         for (j = 0; j < 2; ++j) {
            if (mp->bond_array[i].atoms[j] < 1 || mp->bond_array[i].atoms[j] > mp->n_atoms)
            {
               snprintf(msg_buffer, MAXMSG,
                  "%10s    : illegal atom # (%d, max allowed is %d) in bond %d",
                  mp->name, mp->bond_array[i].atoms[j], mp->n_atoms, i + 1);
               AddMsgToList(msg_buffer);
               result |= SIZE_CHECK_FAILED;
            }
         }
      }
   }

   if ((result & SIZE_CHECK_FAILED) == 0)
   {
      if (convert_atom_texts)
      {
         tmp = ConvertAtomAliases(mp);
         if (tmp == 0) result |= ALIAS_CONVERSION_FAILED;
         if (tmp == 1) result |= TRANSFORMED;
      }

      if (convert_stext)
	 new_data_list = ConvertSTEXTToData(mp, new_data_list);

      if (trans_pairs)
      {
	 tmp = TransformAugmentedAtoms(mp,trans_pairs,ntrans);
	 if (tmp) result |= TRANSFORMED;
      }

      stereo_result = DubiousStereochemistry(mp);
      if (FixDubious3DMolecule(mp) & CONVERTED_TO_2D)
      {
         stereo_result = TRUE;
	 result |= DUBIOUS_STEREO_REMOVED;
      }

      if (remove_minor_fragments)
      {
         /* Add MF_PRE and MW_PRE data fields */
         new_data_list = AddMWMF(new_data_list, mp, "MW_PRE", "MF_PRE");
	 new_mp = StripSmallFragments(CopyMolecule(mp), &fragments_found);
         if (new_mp)      // new molecule data structure has been allocated
         {                // => need to overwrite components pointed to by input pointer
                          //    and deallocate skeleton container object
             // deallocate objects pointed to by components of input (*mp)
             FreeMoleculeChildObjects(mp);
             // shallow copy (*new_mp) into 'emptied' container
             (*mp) = (*new_mp);
             // deallocate new container struct
             MyFree((char *) new_mp);
         }
	 if (fragments_found) result |= FRAGMENTS_FOUND;
         /* Add MF_POST and MW_POST data fields */
         new_data_list = AddMWMF(new_data_list, mp, "MW_POST", "MF_POST");
      }

      for (i=0; i<ntautomers; i++)         /* do tautomer standardization */
      {
// fprintf(stderr, "tautomerizing with rule %d\n", i);
         for (j=0;j<3;j++)      // limit to 3 run per rule
         {
            tmp = ApplyTautomer(mp, from_tautomer[i], to_tautomer[i]);
            if (!tmp) break;
            result |= TAUTOMER_TRANSFORMED;
	    snprintf(msg_buffer, MAXMSG,
		    "%10s: has been tautomerized with rule '%s'",
		    mp->name, from_tautomer[i]->name);
	    AddMsgToList(msg_buffer);
         }
      }

      if (!IsNULL(data_list)  &&  !IsNULL(new_data_list))
      {        // append new data list if any
          for (dph=data_list; !IsNULL(dph->next); dph=dph->next)
              ;
          dph->next = new_data_list;
      }

      if (stereo_result == EITHER_BOND_FOUND)   /* looks for EITHER bonds */
      {
	 result |= EITHER_WARNING;
	 if (stereo_file)
	    Cinderella(stereo_file, mp, data_list, result_as_data, "DUBIOUS_STEREO");
	 if (log_file) PrintMsgs(log_file);
	 RemoveDubiousStereochemistry(mp);
	 result |= DUBIOUS_STEREO_REMOVED;
      }
      else if (stereo_result > EITHER_BOND_FOUND) /* more severe errors */
      {
	 result |= STEREO_ERROR;
	 if (check_stereo)
	    result |= BAD_MOLECULE;
	 else
	 {
	    if (stereo_file)
	       Cinderella(stereo_file, mp, data_list, result_as_data, "DUBIOUS_STEREO");
	    if (log_file) PrintMsgs(log_file);
	    RemoveDubiousStereochemistry(mp);
	    result |= DUBIOUS_STEREO_REMOVED;
	 }
      }

      if (check_collision && AtomClash(mp))
	 result |= ATOM_CLASH;
      if (good_atoms && !CheckAtoms(mp, good_atoms, ngood))
	 result |= ATOM_CHECK_FAILED;
      if (check_stereo && !CheckStereo(mp))
	 result |= STEREO_ERROR;

      if (charges_read && TotalCharge(mp) != desired_charge)
      {
// fprintf(stderr,"recharging '%s'\n",mp->name);
	 tmp = RechargeMolecule(mp, desired_charge, &ndeprot, &nrefine);
	 if (mp->symbol_lists || mp->prop_lines) tmp = FALSE; 
	 charge_bad = !tmp;
	 if (!charge_bad) result |= RECHARGED;
	 else             result |= BAD_MOLECULE;
      }

      if (groups_to_sgroups) ConvertGroupsToSGroups(mp);

      stereo_bad = FALSE;
      for (i=0; i<nstereopat; i++)
      {
	 ssp = stereo_patterns[i];
	 tmp = ForceStereoTemplate(mp, ssp);
	 if (tmp == (-1))
	 {
	    result |= STEREO_FORCED_BAD;
	    snprintf(msg_buffer, MAXMSG,
		    "%10s: problem enforcing stereochemistry of '%s'",
		    mp->name, ssp->name);
	    AddMsgToList(msg_buffer);
	 }
	 else if (tmp == 15)
	 {
	    result |= STEREO_TRANSFORMED;
	    snprintf(msg_buffer, MAXMSG,
		    "%10s: stereochemistry of '%s' enforced",
		    mp->name, ssp->name);
	    AddMsgToList(msg_buffer);
	    if (stereo_file)
	       Cinderella(stereo_file, mp, data_list, result_as_data, "STEREO_FORCED");
	 }
      }

      for (i=0; i<npat; i++)         /* do template cleaning */
      {
	 ssp = patterns[i];
// fprintf(stderr, "cleaning with pattern %d\n", i);
	 tmp = TemplateClean(mp, ssp);
	 if (tmp)
	 {
	    result |= TEMPLATE_TRANSFORMED;
	    snprintf(msg_buffer, MAXMSG,
		    "%10s: has been cleaned with template '%s'",
		    mp->name, ssp->name);
	    AddMsgToList(msg_buffer);
	 }
      }

      for (i=0; i<nrpat; i++)         /* do template rotation */
      {
	 ssp = rotate_patterns[i];
// fprintf(stderr, "rotating with pattern %d\n", i);
	 tmp = TemplateRotate(mp, ssp);
	 if (tmp)
	 {
	    result |= TEMPLATE_TRANSFORMED;
	    snprintf(msg_buffer, MAXMSG,
		    "%10s: has been rotated by template '%s'",
		    mp->name, ssp->name);
	    AddMsgToList(msg_buffer);
	 }
      }
   }

   if (log_file) FlushMsgs(log_file);

   if (IsNULL(data_list))  // there was no input list
   {                               // => free the new items
      while (!IsNULL(new_data_list))
      {
         data_list = new_data_list->next;
         MyFree((char *)new_data_list);
         new_data_list = data_list;
      }
   }

   if (log_file)
   {
      FlushMsgs(log_file);
      fflush(log_file);
   }
   return (result);
}

/**
 * Wrapper around RunStruchk() to be called from Oracle call-outs.
 * Note: InitCheckMol() may have been called before to initialize
 * options.
 *
 * Note 2: This architecture does not work if the structure that mol
 * points to is replaced! This may be the case when removing fragments.
 */
int CheckMol(struct reaccs_molecule_t *mol)
{
   int lastDontExit, tmp;
   lastDontExit = dontExit;
   dontExit = TRUE;
   tmp = RunStruchk(&mol, (struct data_line_t *)NULL);
   dontExit = lastDontExit;
   return tmp;
}

/**************************************************************
 * CloseOpenFiles closes the open global files
 **************************************************************/
void CloseOpenFiles(void)
{
   if (bad_file)
   {
      fclose(bad_file);
      bad_file = (FILE *)NULL;
   }
   if (stereo_file)
   {
      fclose(stereo_file);
      stereo_file = (FILE *)NULL;
   }
   if (log_file)
   {
      fclose(log_file);
      log_file = (FILE *)NULL;
   }
   if (aa_log)
   {
      fclose(aa_log);
      aa_log = (FILE *)NULL;
   }
}

#ifdef MAIN

int main(int argc, char *argv[])
{
   int is_mol_file;
   struct reaccs_molecule_t *mp, *old_mp;
   struct data_line_t *data_list, *dph;

   int flags;

   int n_bad_molecules         = 0,
       n_good_molecules        = 0,
       n_transformed_molecules = 0,
       n_dubious               = 0;

   int errflg;

// FORTIFY long bytes_allocated;
// FORTIFY Fortify_SetOutputFunc(logFortifyMessage);
// FORTIFY Fortify_OutputStatistics();
// FORTIFY fprintf(stderr, "fortify enabled\n");
   
#ifdef unix
   profile_name = ".struchkrc";
#else
   profile_name = "struchk.rc";
#endif

   profile = fopen(profile_name, "r");

   errflg = Initialize(profile, argc, argv);

   if (IsNULL(outfile)) outfile = stdout; 

   if (print_versions)
   {
      fprintf(stdout,"STRUCHK version %s\n",struchk_version);
      if (*aa_trans_version)
         fprintf(stdout, "augmented atom transformation version %s\n", aa_trans_version);
      if (*aa_check_version)
         fprintf(stdout, "augmented atom check version %s\n", aa_check_version);
      if (*pKa_version)
         fprintf(stdout, "acidity table version %s\n", pKa_version);
   }

   if (errflg)  /* illegal option -> show usage */
      MessageExit( "usage: %s {options}\n", argv[0], (char *)NULL);

   if (IsNULL(fp))
   {
      fprintf(stderr,"No input file given\n");
      MessageExit( "usage: %s {options}\n", argv[0], (char *)NULL);
   }

// FORTIFY Fortify_EnterScope();
// FORTIFY bytes_allocated = fortifySet();
// FORTIFY fprintf(stderr, "after entering scope\n");
// FORTIFY Fortify_OutputStatistics();

   mp = TypeAlloc(1,struct reaccs_molecule_t);
   while (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,""))
   {
      old_mp = CopyMolecule(mp);
      is_mol_file = (fp->status == FORTRAN_EOF);

      data_list = ReadMACCSDataLines(fp);

      if (!is_mol_file) GetBuffer(fp);                /* Skip '$$$$'-line */

      flags = RunStruchk(&mp, data_list);
      if (flags & DUBIOUS_STEREO_REMOVED) n_dubious++;

      if (key[0])                           /* look for key field */
         for (dph=data_list; !IsNULL(dph); dph=dph->next)
         {
            if (IsFieldHeader(dph->data, key))
            {
               if (dph->next->data[0])
                  strcpy(mp->name, dph->next->data);
               if (squeeze_identifiers) Squeeze(mp->name);
	       dph=dph->next;
	    }
         }

      if (flags & BAD_SET)
      {
         n_bad_molecules++;
         if (bad_file)
	    Cinderella(bad_file, old_mp, data_list,
		       result_as_data, "BAD");
         if (!remove_bad_mols)
            Cinderella(outfile, old_mp, data_list,
                       result_as_data, "BAD");
      }
      else
      {
         if (flags & TRANSFORMED_SET)
         {
            Cinderella(outfile, mp, data_list,
                       result_as_data, "TRANSFORMED");
            n_good_molecules++;
            n_transformed_molecules++;
         }
         else
         {
            Cinderella(outfile, mp, data_list,
                       result_as_data, "GOOD");
            n_good_molecules++;
         }
      }
      FlushMsgs(log_file);

      FreeMolecule(mp);
      FreeMolecule(old_mp);
      while (!IsNULL(data_list))
      {
         dph = data_list->next;
         MyFree((char *)data_list);
         data_list = dph;
      }
      mp = TypeAlloc(1,struct reaccs_molecule_t);
   }

// Code for pKa optimization
// PrintChargeFooter();
//

   fprintf(stderr, " %4d good molecules, %d of which have been transformed\n",
           n_good_molecules, n_transformed_molecules);
   fprintf(stderr," %4d bad molecules\n", n_bad_molecules);
   fprintf(stderr, " %4d molecules with dubious stereochemistry\n", n_dubious);

   FreeMolecule(mp);
   FortranClose(fp); fp = NULL;
   if (outfile != stdout)
   {
      if (!IsNULL(outfile)) fclose(outfile);
   }
   outfile = (FILE *)NULL;
   CloseOpenFiles();

// FORTIFY fprintf(stderr, "before leaving scope\n");
// FORTIFY Fortify_OutputStatistics();
// FORTIFY fortifyTest(bytes_allocated, "Closing struchk");
// FORTIFY fprintf(stderr, "fortify checked\n");
// FORTIFY Fortify_LeaveScope();
   return (EXIT_SUCCESS);
}

#endif
