/*==============================================*/
/* Copyright (c)  2011-2019  NextMove Software  */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#include <stdlib.h>
#include <stdio.h>

#include "toolkit.h"
#include "molhash.h"

#ifdef _WIN32
#define getc_unlocked _fgetc_nolock
#endif

static unsigned int hash_func = 0;
static bool title_only = false;
static bool all_flag = false;
static const char *inpname;
static const char *outname;

#define MAXNORM 10
static unsigned int normidx = 0;
static unsigned int normalizations[MAXNORM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static void ProcessComponent(NMS_pMOL mol, FILE *fp)
{ 
  std::string title;
  std::string hash;
  std::string smi;

  if (!title_only)
    NMS_GENERATE_SMILES(mol, smi);

  hash = MolHash(mol, hash_func);

  fprintf(fp, "%s", hash.c_str());
  if (!title_only)
    fprintf(fp, "\t%s", smi.c_str());
  title = NMS_MOL_GET_TITLE(mol);
  if (!title.empty())
    fprintf(fp, "\t%s", title.c_str());
  fputc('\n', fp);
}


static void ProcessMolecule(NMS_pMOL mol, FILE *fp)
{
  if (all_flag)
    ProcessComponent(mol, fp);
  else {
    std::vector<NMS_MOL*> vmol;
    SplitMolecule(mol, vmol);
    unsigned int largest = 0;
    unsigned int maxsize = 0;
    for(unsigned int i=0; i<vmol.size(); ++i) {
      NMS_pMOL pmol = NMS_MOL_TO_pMOL(*vmol[i]);
      if (NMS_MOL_GET_NUMATOMS(pmol) > maxsize) {
        maxsize = NMS_MOL_GET_NUMATOMS(pmol);
        largest = i;
      }
    }
    ProcessComponent(NMS_MOL_TO_pMOL(*vmol[largest]), fp);
    for (unsigned int i = 0; i < vmol.size(); ++i)
      delete vmol[i];
  }
}


/* Read a line, returning first n characters */
static bool SWReadLine2(FILE *fp, char *buffer, unsigned int n)
{
  char *end = buffer + n;
  char *ptr;
  int ch;

  ptr = buffer;
  do {
    ch = getc_unlocked(fp);
    if (ch == '\n') {
      *ptr = '\0';
      return true;
    }
    if (ch == '\r') {
      *ptr = '\0';
      ch = getc_unlocked(fp);
      if (ch != '\n') {
        if (ch == -1)
          return false;
        ungetc(ch, fp);
      }
      return true;
    }
    if (ch == -1) {
      *ptr = '\0';
      return false;
    }
    *ptr++ = ch;
  } while (ptr < end);
  *ptr = 0;

  /* skip to the end of the line! */
  for (;;) {
    ch = getc_unlocked(fp);
    if (ch == '\n')
      return true;
    if (ch == '\r') {
      ch = getc_unlocked(fp);
      if (ch != '\n') {
        if (ch == -1)
          return false;
        ungetc(ch, fp);
      }
      return true;
    }
    if (ch == -1)
      return false;
  }
}

static void ProcessFile(FILE *ifp, FILE *ofp)
{
  char buffer[8192];
  while (SWReadLine2(ifp, buffer, 8190)) {
    NMS_MOL* mol = NMS_SMILES_TO_MOL(buffer);
    if (mol != (NMS_MOL*)0) {
      NMS_pMOL pmol = NMS_MOL_TO_pMOL(*mol);
      for (unsigned int i=0; i<normidx; ++i)
        Strip(pmol, normalizations[i]);
      ProcessMolecule(pmol, ofp);
      delete mol;
    }
  }
}

unsigned int HashOption(const char *ptr)
{
  if (!strcmp(ptr, "-g")) return HashFunction::AnonymousGraph;
  if (!strcmp(ptr, "-e")) return HashFunction::ElementGraph;
  if (!strcmp(ptr, "-s")) return HashFunction::CanonicalSmiles;
  if (!strcmp(ptr, "-m")) return HashFunction::MurckoScaffold;
  if (!strcmp(ptr, "-em")) return HashFunction::ExtendedMurcko;
  if (!strcmp(ptr, "-mf")) return HashFunction::MolFormula;
  if (!strcmp(ptr, "-ab")) return HashFunction::AtomBondCounts;
  if (!strcmp(ptr, "-dv")) return HashFunction::DegreeVector;
  if (!strcmp(ptr, "-me")) return HashFunction::Mesomer;
  if (!strcmp(ptr, "-ht")) return HashFunction::HetAtomTautomer;
  if (!strcmp(ptr, "-hp")) return HashFunction::HetAtomProtomer;
  if (!strcmp(ptr, "-rp")) return HashFunction::RedoxPair;
  if (!strcmp(ptr, "-ri")) return HashFunction::Regioisomer;
  if (!strcmp(ptr, "-nq")) return HashFunction::NetCharge;
  if (!strcmp(ptr, "-br")) return HashFunction::SmallWorldIndexBR;
  if (!strcmp(ptr, "-brl")) return HashFunction::SmallWorldIndexBRL;
  if (!strcmp(ptr, "-sub")) return HashFunction::ArthorSubstructureOrder;
  return 0;
}


void DisplayHashOptions()
{
  fprintf(stderr, "    -g   anonymous graph [default]\n");
  fprintf(stderr, "    -e   element graph\n");
  fprintf(stderr, "    -s   canonical smiles\n");
  fprintf(stderr, "    -m   Murcko scaffold\n");
  fprintf(stderr, "    -em  Extended Murcko scaffold\n");
  fprintf(stderr, "    -mf  molecular formula\n");
  fprintf(stderr, "    -ab  atom and bond counts\n");
  fprintf(stderr, "    -dv  degree vector\n");
  fprintf(stderr, "    -me  mesomer\n");
  fprintf(stderr, "    -ht  hetatom tautomer\n");
  fprintf(stderr, "    -hp  hetatom protomer\n");
  fprintf(stderr, "    -rp  redox-pair\n");
  fprintf(stderr, "    -ri  regioisomer\n");
  fprintf(stderr, "    -nq  net charge\n");
  fprintf(stderr, "    -sub Arthor substructure order\n");
}

static void DisplayUsage()
{
  fputs("usage:  molhash [options] <infile> [<outfile>]\n",stderr);
  fputs("    Use a hyphen for <infile> to read from stdin\n",stderr);
  fputs("options:\n", stderr);
  fputs("    -a  Process all the molecule (and not just the single largest component)\n", stderr);
  fputs("    -sa Suppress atom stereo\n", stderr);
  fputs("    -sb Suppress bond stereo\n", stderr);
  fputs("    -sh Suppress explicit hydrogens\n", stderr);
  fputs("    -si Suppress isotopes\n", stderr);
  fputs("    -sm Suppress atom maps\n", stderr);
  fputs("    -t  Store titles only\n",stderr);
  fputs("hash type:\n", stderr);
  DisplayHashOptions();
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr;
  int i,j;

  inpname = (const char*)0;
  outname = (const char*)0;
  title_only = false;
  all_flag = false;
  hash_func = 0;

  j = 0;
  for (i=1; i<argc; i++) {
    ptr = argv[i];
    if (ptr[0]=='-' && ptr[1]) {
      if (!strcmp(ptr, "-a")) {
        all_flag = true;
      } else if (normidx<MAXNORM && !strcmp(ptr, "-sh")) {
        normalizations[normidx++] = StripType::Hydrogen;
      }
      else if (normidx < MAXNORM && !strcmp(ptr, "-si")) {
        normalizations[normidx++] = StripType::Isotope;
      }
      else if (normidx < MAXNORM && !strcmp(ptr, "-sa")) {
        normalizations[normidx++] = StripType::AtomStereo;
      }
      else if (normidx < MAXNORM && !strcmp(ptr, "-sb")) {
        normalizations[normidx++] = StripType::BondStereo;
      }
      else if (normidx < MAXNORM && !strcmp(ptr, "-sm")) {
        normalizations[normidx++] = StripType::AtomMap;
      } else if (!strcmp(ptr,"-t")) {
        title_only = true;
      } else {
        hash_func = HashOption(ptr);
        if (hash_func == 0)
          DisplayUsage();
      }
    } else switch (j++) {
    case 0:  inpname = ptr;  break;
    case 1:  outname = ptr;  break;
    case 2:  DisplayUsage();
    }
  }

  if (j < 1)
    DisplayUsage();
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc,argv);

  FILE *ifp;
  if (inpname && *inpname && strcmp(inpname,"-")) {
    ifp = fopen(inpname,"rb");
    if (!ifp) {
      fprintf(stderr,"Error: Cannot read input file: %s\n",inpname);
      exit(1);
    }
  } else ifp = stdin;

  FILE *ofp;
  if (outname && *outname && strcmp(outname,"-")) {
    ofp = fopen(outname,"wb");
    if (!ofp) {
      fprintf(stderr,"Error: Cannot create output file: %s\n",outname);
      exit(1);
    }
  } else ofp = stdout;

  ProcessFile(ifp,ofp);

  if (ofp != stdout)
    fclose(ofp);
  if (ifp != stdin)
    fclose(ifp);
  return 0;
}

