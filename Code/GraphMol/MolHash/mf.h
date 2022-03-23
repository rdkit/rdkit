//
//  Copyright (C) 2016-2022 NextMove Software and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
#ifndef NMS_MOLFORMULA_H
#define NMS_MOLFORMULA_H

static unsigned char OrganicHillOrder[119] = {
    6,   /* C  */
    1,   /* H  */
    89,  /* Ac */
    47,  /* Ag */
    13,  /* Al */
    95,  /* Am */
    18,  /* Ar */
    33,  /* As */
    85,  /* At */
    79,  /* Au */
    5,   /* B  */
    56,  /* Ba */
    4,   /* Be */
    107, /* Bh */
    83,  /* Bi */
    97,  /* Bk */
    35,  /* Br */
    20,  /* Ca */
    48,  /* Cd */
    58,  /* Ce */
    98,  /* Cf */
    17,  /* Cl */
    96,  /* Cm */
    112, /* Cn */
    27,  /* Co */
    24,  /* Cr */
    55,  /* Cs */
    29,  /* Cu */
    105, /* Db */
    110, /* Ds */
    66,  /* Dy */
    68,  /* Er */
    99,  /* Es */
    63,  /* Eu */
    9,   /* F  */
    26,  /* Fe */
    114, /* Fl */
    100, /* Fm */
    87,  /* Fr */
    31,  /* Ga */
    64,  /* Gd */
    32,  /* Ge */
    2,   /* He */
    72,  /* Hf */
    80,  /* Hg */
    67,  /* Ho */
    108, /* Hs */
    53,  /* I  */
    49,  /* In */
    77,  /* Ir */
    19,  /* K  */
    36,  /* Kr */
    57,  /* La */
    3,   /* Li */
    103, /* Lr */
    71,  /* Lu */
    116, /* Lv */
    115, /* Mc */
    101, /* Md */
    12,  /* Mg */
    25,  /* Mn */
    42,  /* Mo */
    109, /* Mt */
    7,   /* N  */
    11,  /* Na */
    41,  /* Nb */
    60,  /* Nd */
    10,  /* Ne */
    113, /* Nh */
    28,  /* Ni */
    102, /* No */
    93,  /* Np */
    8,   /* O  */
    118, /* Og */
    76,  /* Os */
    15,  /* P  */
    91,  /* Pa */
    82,  /* Pb */
    46,  /* Pd */
    61,  /* Pm */
    84,  /* Po */
    59,  /* Pr */
    78,  /* Pt */
    94,  /* Pu */
    88,  /* Ra */
    37,  /* Rb */
    75,  /* Re */
    104, /* Rf */
    111, /* Rg */
    45,  /* Rh */
    86,  /* Rn */
    44,  /* Ru */
    16,  /* S  */
    51,  /* Sb */
    21,  /* Sc */
    34,  /* Se */
    106, /* Sg */
    14,  /* Si */
    62,  /* Sm */
    50,  /* Sn */
    38,  /* Sr */
    73,  /* Ta */
    65,  /* Tb */
    43,  /* Tc */
    52,  /* Te */
    90,  /* Th */
    22,  /* Ti */
    81,  /* Tl */
    69,  /* Tm */
    117, /* Ts */
    92,  /* U  */
    23,  /* V  */
    74,  /* W  */
    0,   /* X  */
    54,  /* Xe */
    39,  /* Y  */
    70,  /* Yb */
    30,  /* Zn */
    40   /* Zr */
};

static unsigned char InorganicHillOrder[119] = {
    89,  /* Ac */
    47,  /* Ag */
    13,  /* Al */
    95,  /* Am */
    18,  /* Ar */
    33,  /* As */
    85,  /* At */
    79,  /* Au */
    5,   /* B  */
    56,  /* Ba */
    4,   /* Be */
    107, /* Bh */
    83,  /* Bi */
    97,  /* Bk */
    35,  /* Br */
    6,   /* C  */
    20,  /* Ca */
    48,  /* Cd */
    58,  /* Ce */
    98,  /* Cf */
    17,  /* Cl */
    96,  /* Cm */
    112, /* Cn */
    27,  /* Co */
    24,  /* Cr */
    55,  /* Cs */
    29,  /* Cu */
    105, /* Db */
    110, /* Ds */
    66,  /* Dy */
    68,  /* Er */
    99,  /* Es */
    63,  /* Eu */
    9,   /* F  */
    26,  /* Fe */
    114, /* Fl */
    100, /* Fm */
    87,  /* Fr */
    31,  /* Ga */
    64,  /* Gd */
    32,  /* Ge */
    1,   /* H  */
    2,   /* He */
    72,  /* Hf */
    80,  /* Hg */
    67,  /* Ho */
    108, /* Hs */
    53,  /* I  */
    49,  /* In */
    77,  /* Ir */
    19,  /* K  */
    36,  /* Kr */
    57,  /* La */
    3,   /* Li */
    103, /* Lr */
    71,  /* Lu */
    116, /* Lv */
    115, /* Mc */
    101, /* Md */
    12,  /* Mg */
    25,  /* Mn */
    42,  /* Mo */
    109, /* Mt */
    7,   /* N  */
    11,  /* Na */
    41,  /* Nb */
    60,  /* Nd */
    10,  /* Ne */
    113, /* Nh */
    28,  /* Ni */
    102, /* No */
    93,  /* Np */
    8,   /* O  */
    118, /* Og */
    76,  /* Os */
    15,  /* P  */
    91,  /* Pa */
    82,  /* Pb */
    46,  /* Pd */
    61,  /* Pm */
    84,  /* Po */
    59,  /* Pr */
    78,  /* Pt */
    94,  /* Pu */
    88,  /* Ra */
    37,  /* Rb */
    75,  /* Re */
    104, /* Rf */
    111, /* Rg */
    45,  /* Rh */
    86,  /* Rn */
    44,  /* Ru */
    16,  /* S  */
    51,  /* Sb */
    21,  /* Sc */
    34,  /* Se */
    106, /* Sg */
    14,  /* Si */
    62,  /* Sm */
    50,  /* Sn */
    38,  /* Sr */
    73,  /* Ta */
    65,  /* Tb */
    43,  /* Tc */
    52,  /* Te */
    90,  /* Th */
    22,  /* Ti */
    81,  /* Tl */
    69,  /* Tm */
    117, /* Ts */
    92,  /* U  */
    23,  /* V  */
    74,  /* W  */
    0,   /* X  */
    54,  /* Xe */
    39,  /* Y  */
    70,  /* Yb */
    30,  /* Zn */
    40   /* Zr */
};

//  x0    x1    x2    x3    x4    x5    x6    x7    x8    x9
static const char *symbol[119] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",   //   x
    "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",   //  1x
    "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",  //  2x
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",   //  3x
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",  //  4x
    "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",  //  5x
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",  //  6x
    "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  //  7x
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",  //  8x
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",  //  9x
    "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",  // 10x
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

#endif  // NMS_MOLFORMULA_H
