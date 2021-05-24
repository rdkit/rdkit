//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#if defined(__CYGWIN__) && !defined(_GNU_SOURCE)
// -std=c++11 doesn't declare strtok_r
#define _GNU_SOURCE
#endif

#include <cstring>
#include <cctype>
#include "Pattern.h"

namespace RDKit {
namespace StructureCheck {

static const char *AtomSymbol[] = {
    // Periodic Table
    "*",   // ANY_ELEMENT = 0,
    "H",   // 1
    "He",  // 2
    "Li", "Be", "B",
    "C",  // 6
    "N",  "O",  "F",
    "Ne",  // 10
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl",
    "Ar",  // 18
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
    "Ni",  // 28
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",  // 36
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd",  // 46
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
    "Xe",  // 54
    "Cs", "Ba",

    "La",  // 57
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb",
    "Lu",  // 71

    "Hf", "Ta", "W",  "Re", "Os", "Ir",
    "Pt",  // 78
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn",  // 86
    "Fr", "Ra",

    "Ac",  // 89
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No",
    "Lr",  // 103

    "Rf", "Db", "Sg", "Bh", "Hn",
    "Mt"  // 109
};

class AtomSymbolMapper {
  std::map<std::string, unsigned> SymbolMap;

 public:
  AtomSymbolMapper() {
    for (unsigned n = 0; n < 110; n++) SymbolMap[AtomSymbol[n]] = n;
  }
  inline unsigned getAtomicNumber(const std::string symbol) const {
    return SymbolMap.find(symbol)->second;
  }
};
static const AtomSymbolMapper smap;

unsigned getAtomicNumber(const std::string symbol) {
  return smap.getAtomicNumber(symbol);
}

// predefined generic atom type sets for use in STRUCHK
static const char *HC_table[] = /* pseudosymbol "G" */
    {"H", "C", nullptr};

static const char *non_metal_hetero_elements[] = /* pseudosymbol "Q" */
    {
        "He",    "B",  "N",  "O",  "F",  "Ne", "Si", "P",  "S",  "Cl", "Ar",
        "As",    "Se", "Br", "Kr", "Sb", "Te", "I",  "Xe", "At", /* "Rn", This
                                                                    element must
                                                                    be
                                                                    removed */
        nullptr, /* because of a trick in utils.c */
};

static const char *metals[] = /* pseudosymbol "M" */
    {
        "Li", "Be", "Na", "Mg", "Al",    "K",  "Ca", "Sc", "Ti", "V",  "Cr",
        "Mn", "Fe", "Co", "Ni", "Cu",    "Zn", "Ga", "Rb", "Sr", "Y",  "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh",    "Pd", "Ag", "Cd", "In", "Sn", "Cs",
        "Ba", "La", "Ce", "Pr", "Nd",    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu",    "Hf", "Ta", "W",  "Re", "Os", "Ir",
        "Pt", "Au", "Hg", "Tl", "Pb",    "Bi", "Po", "Fr", "Ra", "Ac", "Th",
        "Pa", "U",  "Np", "Pu", nullptr,
};

static const char *non_metal_small_solution[] = /* pseudosymbol "Qs" */
    {
        "H", "B", "C",  "N",  "O",  "F", "Si",
        "P", "S", "Cl", "Se", "Br", "I", nullptr,
};

static const char *alkali_metals[] = /* pseudosymbol "alk" */
    {
        "Li", "Na", "K", "Rb", "Cs", "Fr", nullptr,
};

static const char *gr2[] = /* pseudosymbol "gr2" */
    {
        "Be", "Mg", "Ca", "Sr", "Ba", "Ra", nullptr,
};

static const char *gr3[] = /* pseudosymbol "gr3" */
    {
        "B", "Al", "Ga", "In", "Tl", nullptr,
};

static const char *gr4[] = /* pseudosymbol "gr4" */
    {
        "C", "Si", "Ge", "Sn", "Pb", nullptr,
};

static const char *ONS_table[] = /* pseudosymbol "ONS" or "ons" */
    {"O", "N", "S", nullptr};

static const char *on2[] = /* pseudosymbol "on2" */
    {
        "O", "N", "S", "P", "Se", "Te", "Po", nullptr,
};

static const char *halogenes[] = /* pseudosymbol "X" or "hal" */
    {"F", "Cl", "Br", "I", "At", nullptr};

static const char *ha2[] = /* pseudosymbol "ha2" */
    {"Cl", "Br", "I", "At", nullptr};

static const char *transition_metals[] = /* pseudosymbol "trn" */
    {
        "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",    "Zn", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",    "La", "Hf",
        "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", nullptr,
};

static const char *tra[] = /* pseudosymbol "tra" */
    {
        "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Zr", "Nb", "Mo", "Tc",
        "Ru", "Rh", "Pd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", nullptr,
};

static const char *trb[] = /* pseudosymbol "trb" */
    {
        "Cu", "Zn", "Ag", "Cd", "Au", "Hg", nullptr,
};

static const char *tm1[] = /* pseudosymbol "tm1" */
    {
        "Cu",
        "Ag",
        "Au",
        nullptr,
};

static const char *tm2[] = /* pseudosymbol "tm2" */
    {
        "Zn",
        "Cd",
        "Hg",
        nullptr,
};

static const char *tm3[] = /* pseudosymbol "tm3" */
    {
        "Sc",
        "Y",
        "La",
        nullptr,
};

static const char *tm4[] = /* pseudosymbol "tm4" */
    {
        "Ti",
        "Zr",
        "Hf",
        nullptr,
};

static const char *tm5[] = /* pseudosymbol "tm5" */
    {
        "V",
        "Nb",
        "Ta",
        nullptr,
};

static const char *tm6[] = /* pseudosymbol "tm6" */
    {
        "Cr",
        "Mo",
        "W",
        nullptr,
};

static const char *tm7[] = /* pseudosymbol "tm7" */
    {
        "Mn",
        "Tc",
        "Re",
        nullptr,
};

static const char *tm8[] = /* pseudosymbol "tm8" */
    {
        "Fe", "Co", "Ni", "Ru", "Rh", "Pd", "Os", "Ir", "Pt", nullptr,
};

static const char *lanthanoids[] = /* pseudosymbol "lan" */
    {
        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",    "Tb",
        "Dy", "Ho", "Er", "Tm", "Yb", "Lu", nullptr,
};

static const char *amino_acids[] = /* pseudosymbol "Ami" or "ami"*/
    {
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu",
        "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
        "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", nullptr,
};

static bool IsInStringTable(const char *symbol, const char *table[]) {
  // Checks if the string symbol is listed in table[]
  const char **stringp;
  for (stringp = table; *stringp; stringp++)
    if (0 == strcmp(*stringp, symbol)) return true;
  return false;
}

bool AtomSymbolMatch(const std::string symbol, const std::string pattern) {
  /*
   * Returns TRUE if symbol is in the comma separated list of atom symbols
   * stored in pattern and FALSE otherwise.
   * There are also a number of standard atom type lists like "alk" for alkali
   * metals or
   * "Q" for non-C/non-H defined above as arrays of strings.
   */
  char *context;
#ifdef WIN32
#define strtok_r strtok_s  // thread safe strtok()
#endif
  const char *atsym = symbol.c_str();
  char pat_buf[512];
  char *tokp;

  strcpy(pat_buf, pattern.c_str());
  for (tokp = strtok_r(pat_buf, ",", &context); tokp;
       tokp = strtok_r((char *)nullptr, ",", &context)) {
    if (islower(*tokp)) {
      if (0 == strcmp("alk", tokp)) {
        if (IsInStringTable(atsym, alkali_metals)) return true;
      } else if (0 == strcmp("gr2", tokp)) {
        if (IsInStringTable(atsym, gr2)) return true;
      } else if (0 == strcmp("gr3", tokp)) {
        if (IsInStringTable(atsym, gr3)) return true;
      } else if (0 == strcmp("gr4", tokp)) {
        if (IsInStringTable(atsym, gr4)) return true;
      } else if (0 == strcmp("ons", tokp)) {
        if (IsInStringTable(atsym, ONS_table)) return true;
      } else if (0 == strcmp("on2", tokp)) {
        if (IsInStringTable(atsym, on2)) return true;
      } else if (0 == strcmp("hal", tokp)) {
        if (IsInStringTable(atsym, halogenes)) return true;
      } else if (0 == strcmp("ha2", tokp)) {
        if (IsInStringTable(atsym, ha2)) return true;
      } else if (0 == strcmp("trn", tokp)) {
        if (IsInStringTable(atsym, transition_metals)) return true;
      } else if (0 == strcmp("tra", tokp)) {
        if (IsInStringTable(atsym, tra)) return true;
      } else if (0 == strcmp("trb", tokp)) {
        if (IsInStringTable(atsym, trb)) return true;
      } else if (0 == strcmp("tm1", tokp)) {
        if (IsInStringTable(atsym, tm1)) return true;
      } else if (0 == strcmp("tm2", tokp)) {
        if (IsInStringTable(atsym, tm2)) return true;
      } else if (0 == strcmp("tm3", tokp)) {
        if (IsInStringTable(atsym, tm3)) return true;
      } else if (0 == strcmp("tm4", tokp)) {
        if (IsInStringTable(atsym, tm4)) return true;
      } else if (0 == strcmp("tm5", tokp)) {
        if (IsInStringTable(atsym, tm5)) return true;
      } else if (0 == strcmp("tm6", tokp)) {
        if (IsInStringTable(atsym, tm6)) return true;
      } else if (0 == strcmp("tm7", tokp)) {
        if (IsInStringTable(atsym, tm7)) return true;
      } else if (0 == strcmp("tm8", tokp)) {
        if (IsInStringTable(atsym, tm8)) return true;
      } else if (0 == strcmp("lan", tokp)) {
        if (IsInStringTable(atsym, lanthanoids)) return true;
      } else if (0 == strcmp("ami", tokp)) {
        if (IsInStringTable(atsym, amino_acids)) return true;
      }
    }
    if (0 == strcmp(atsym, tokp)) return true;
  }

  if (0 == strcmp("A", pattern.c_str()))
    return (0 != strcmp("H", atsym));
  else if (0 == strcmp("Qs", pattern.c_str()))
    return (IsInStringTable(atsym, non_metal_small_solution));
  else if (0 == strcmp("G", pattern.c_str()))
    return (IsInStringTable(atsym, HC_table));
  else if (0 == strcmp("ONS", pattern.c_str()))
    return (IsInStringTable(atsym, ONS_table));
  else if (0 == strcmp("X", pattern.c_str()))
    return (IsInStringTable(atsym, halogenes));
  else if (0 == strcmp("Q", pattern.c_str()))
    return (IsInStringTable(atsym, non_metal_hetero_elements));
  else if (0 == strcmp("M", pattern.c_str()))
    return (IsInStringTable(atsym, metals));
  else if (0 == strcmp("Ami", pattern.c_str()))
    return (IsInStringTable(atsym, amino_acids));
  return false;
}

}  // namespace StructureCheck
}  // namespace RDKit
