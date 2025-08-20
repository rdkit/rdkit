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
package novartis.chemistry.molecule;

import java.util.StringTokenizer;

/**
 * This class implements static functions to lookup properties
 * in a periodic table. The class name should have been PTableEntry
 * but, for the static functions, PTable.f() looks more natural.
 */
public class PTable
{
   // instance variables
   int    atomic_number;
   String symbol;
   int    normal_valence;
   float  natural_mass;
   float  electronegativity;

   /* list of allowed valences of this element, null means any */
   int allowed_valences[] = null;

   // special atom type constants, i.e pseudo atomic numbers
   final static int HALOGENE     = 120;
   final static int HETERO       = 121;
   final static int METAL        = 122;
   final static int RGROUP       = 123;
   final static int ERROR_NUMBER = (-1);
   final static int QUERY_STRING = (-2);

   /**
    * Constructor periodic table row.
    */
   PTable(String symbol, int valence, double mass, double elneg, int ord)
   {
      this.symbol            = symbol;
      this.normal_valence    = valence;
      this.natural_mass      = (float)mass;
      this.electronegativity = (float)elneg;
      this.atomic_number     = ord;
      this.allowed_valences  = null;
   }

   /**
    * Constructor periodic table row with constrained valence.
    */
   PTable(String symbol, int valence, double mass, double elneg, int ord,
          int from_valence, int to_valence, int step_valence)
   {
      this.symbol            = symbol;
      this.normal_valence    = valence;
      this.natural_mass      = (float)mass;
      this.electronegativity = (float)elneg;
      this.atomic_number     = ord;

      // initialize valence table
      int count = 0;
      for (int i=from_valence; i<=to_valence; i+=step_valence)
         count++;
      this.allowed_valences   = new int[count];
      for (int i=from_valence, j=0; i<=to_valence; i+=step_valence, j++)
         this.allowed_valences[j] = i;
   }

   /**
    * Conversion of atomic symbol to the atomic ordinal number in
    * the periodic table of elements.
    */
   public static int SymbolToAtomicNumber(String symbol)
   {
      for (int i=0; i<ptable.length; i++)
         if (symbol.equals(ptable[i].symbol))
	    return (ptable[i].atomic_number);

      return (ERROR_NUMBER);
   }

   /**
    * Conversion of atomic ordinal number to the element symbol.
    */
   public static String AtomicNumberToSymbol(int atomic_number)
   {
      for (int i=0; i<ptable.length; i++)
         if (atomic_number == ptable[i].atomic_number)
	    return (ptable[i].symbol);

      return (null);
   }

   /**
    * Convert mass difference values to actual isotopic values for a given
    * element defined by its atom symbol.
    */
   public static int mdiffToIsotope(String symbol, int mdiff)
   {  
      if (mdiff == 0) return (0);
      for (int i=0; i<ptable.length; i++)
         if (symbol.equals(ptable[i].symbol))
         {
            return ((int)(ptable[i].natural_mass + mdiff + 0.49999));
         }
      return (0);
   }


   /**
     * Convert element symbol to its
     * standard valence value
     */
   public static int symbolToValence(String symbol)
   {
      for (int i=0; i<ptable.length; i++)
         if (symbol.equals(ptable[i].symbol))
            return (ptable[i].normal_valence);
      return (0);
   }

   // actual static table
   static PTable[] ptable =
      {
         new PTable("C" ,  4,   12.01115,  2.55,   6,  4, 4, 1),
         new PTable("N" ,  3,   14.00670,  3.04,   7,  3, 5, 2),
         new PTable("O" ,  2,   15.99940,  3.44,   8,  2, 2, 1),
         new PTable("H" ,  1,    1.00797,  2.20,   1, -1, 1, 2),
         new PTable("D" ,  1,    2.01400,  2.20,   1),
         new PTable("T" ,  1,    3.01610,  2.20,   1),
         new PTable("He",  0,    4.00260,  0.00,   2),
         new PTable("Li",  1,    6.93900,  0.98,   3, -1, 1, 2),
         new PTable("Be",  2,    9.01220,  1.57,   4, -2, 2, 2),
         new PTable("B" ,  3,   10.81100,  2.04,   5,  3, 3, 1),
         new PTable("C" ,  4,   12.01115,  2.55,   6,  4, 4, 1),
         new PTable("N" ,  3,   14.00670,  3.04,   7,  3, 5, 2),
         new PTable("O" ,  2,   15.99940,  3.44,   8,  2, 2, 1),
         new PTable("F" ,  1,   18.99840,  3.98,   9,  1, 1, 1),
         new PTable("Ne",  0,   20.18300,  0.00,  10),
         new PTable("Na",  1,   22.98980,  0.93,  11, -1, 1, 2),
         new PTable("Mg",  2,   24.31200,  1.31,  12, -2, 2, 2),
         new PTable("Al",  3,   26.98150,  1.61,  13, -3, 3, 2),
         new PTable("Si",  4,   28.08600,  1.90,  14,  4, 4, 1),
         new PTable("P" ,  5,   30.97380,  2.19,  15,  3, 5, 2),
         new PTable("S" ,  2,   32.06400,  2.58,  16,  2, 6, 2),
         new PTable("Cl",  1,   35.45300,  3.16,  17,  1, 7, 2),
         new PTable("Ar",  0,   39.94800,  0.00,  18),
         new PTable("K" ,  1,   39.10200,  0.82,  19, -1, 1, 2),
         new PTable("Ca",  2,   40.08000,  1.00,  20),
         new PTable("Sc",  3,   44.95600,  1.36,  21),
         new PTable("Ti",  3,   47.90000,  1.54,  22),
         new PTable("V" ,  3,   50.94200,  1.63,  23),
         new PTable("Cr",  3,   51.99600,  1.66,  24),
         new PTable("Mn",  4,   54.93800,  1.55,  25),
         new PTable("Fe",  2,   55.84700,  1.83,  26),
         new PTable("Co",  2,   58.93320,  1.88,  27),
         new PTable("Ni",  2,   58.71000,  1.91,  28),
         new PTable("Cu",  1,   63.54600,  1.90,  29),
         new PTable("Zn",  2,   65.37000,  1.65,  30),
         new PTable("Ga",  2,   69.72000,  1.81,  31, -3, 3, 2),
         new PTable("Ge",  4,   72.59000,  2.01,  32),
         new PTable("As",  3,   74.92160,  2.18,  33),
         new PTable("Se",  4,   78.96000,  2.55,  34,  2, 6, 2),
         new PTable("Br",  1,   79.90400,  2.96,  35,  1, 7, 2),
         new PTable("Kr",  0,   83.80000,  0.00,  36),
         new PTable("Rb",  1,   85.47000,  0.82,  37, -1, 1, 2),
         new PTable("Sr",  2,   87.62000,  0.95,  38),
         new PTable("Y" ,  3,   88.90500,  1.22,  39),
         new PTable("Zr",  4,   91.22000,  1.33,  40),
         new PTable("Nb",  3,   92.90600,  1.60,  41),
         new PTable("Mo",  4,   95.94000,  2.16,  42),
         new PTable("Tc",  6,   98.90620,  1.90,  43),
         new PTable("Ru",  4,  101.07000,  2.20,  44),
         new PTable("Rh",  3,  102.90500,  2.28,  45),
         new PTable("Pd",  4,  106.40000,  2.20,  46),
         new PTable("Ag",  1,  107.86800,  1.93,  47),
         new PTable("Cd",  2,  112.40000,  1.69,  48),
         new PTable("In",  3,  114.82000,  1.78,  49, -3, 3, 2),
         new PTable("Sn",  2,  118.69000,  1.96,  50),
         new PTable("Sb",  3,  121.75000,  2.05,  51),
         new PTable("Te",  4,  127.60000,  2.10,  52,  2, 6, 2),
         new PTable("I" ,  1,  126.90440,  2.66,  53,  1, 7, 2),
         new PTable("Xe",  0,  131.30000,  0.00,  54),
         new PTable("Cs",  1,  132.90500,  0.79,  55, -1, 1, 2),
         new PTable("Ba",  2,  137.33000,  0.89,  56),
         new PTable("La",  3,  138.91000,  1.10,  57, -3, 3, 2),
         new PTable("Ce",  3,  140.12000,  1.12,  58, -3, 3, 2),
         new PTable("Pr",  3,  140.90700,  1.13,  59, -3, 3, 2),
         new PTable("Nd",  3,  144.24000,  1.14,  60, -3, 3, 2),
         new PTable("Pm",  3,  145.00000,  1.20,  61, -3, 3, 2),
         new PTable("Sm",  2,  150.35000,  1.17,  62, -3, 3, 2),
         new PTable("Eu",  2,  151.96000,  1.20,  63, -3, 3, 2),
         new PTable("Gd",  3,  157.25000,  1.20,  64, -3, 3, 2),
         new PTable("Tb",  3,  158.92400,  1.20,  65, -3, 3, 2),
         new PTable("Dy",  3,  162.50000,  1.22,  66, -3, 3, 2),
         new PTable("Ho",  3,  164.93000,  1.23,  67, -3, 3, 2),
         new PTable("Er",  3,  167.26000,  1.24,  68, -3, 3, 2),
         new PTable("Tm",  3,  168.93400,  1.25,  69, -3, 3, 2),
         new PTable("Yb",  2,  173.04000,  1.10,  70, -3, 3, 2),
         new PTable("Lu",  3,  174.97000,  1.27,  71, -3, 3, 2),
         new PTable("Hf",  4,  178.49000,  1.30,  72),
         new PTable("Ta",  5,  180.94800,  1.50,  73),
         new PTable("W" ,  6,  183.85000,  2.36,  74),
         new PTable("Re",  0,  186.20000,  1.90,  75),
         new PTable("Os",  3,  190.20000,  2.20,  76),
         new PTable("Ir",  3,  192.20000,  2.20,  77),
         new PTable("Pt",  2,  195.09000,  2.28,  78),
         new PTable("Au",  1,  196.96700,  2.54,  79),
         new PTable("Hg",  2,  200.59000,  2.00,  80),
         new PTable("Tl",  1,  204.37000,  2.04,  81, -3, 3, 2),
         new PTable("Pb",  2,  207.19000,  2.33,  82),
         new PTable("Bi",  3,  208.98000,  2.02,  83),
         new PTable("Po",  0,  209.00000,  2.00,  84),
         new PTable("At",  3,  210.00000,  2.20,  85),
         new PTable("Rn",  0,  222.00000,  0.00,  86),
         new PTable("Fr",  1,  223.00000,  0.70,  87),
         new PTable("Ra",  2,  226.03000,  0.90,  88),
         new PTable("Ac",  0,  227.00000,  1.10,  89),
         new PTable("Th",  4,  232.03800,  1.30,  90),
         new PTable("Pa",  0,  231.04000,  1.50,  91),
         new PTable("U" ,  6,  238.03000,  1.38,  92),
         new PTable("Np",  5,  237.05000,  1.36,  93),
         new PTable("Pu",  4,  244.00000,  1.28,  94),
         new PTable("Am",  4,  243.00000,  1.30,  95),
         new PTable("Cm",  3,  247.00000,  1.30,  96),
         new PTable("Bk",  3,  247.00000,  1.30,  97),
         new PTable("Cf",  0,  251.00000,  1.30,  98),
         new PTable("Es",  0,  254.00000,  1.30,  99),
         new PTable("Fm",  0,  257.00000,  1.30, 100),
         new PTable("Md",  0,  258.00000,  0.00, 101),
         new PTable("No",  0,  259.00000,  0.00, 102),
         new PTable("Lr",  0,  260.00000,  0.00, 103),
         new PTable("X",   1,  260.00000,  0.00, HALOGENE),
         new PTable("Q",   0,  260.00000,  0.00, HETERO),
         new PTable("M",   0,  260.00000,  0.00, METAL),
         new PTable("R#",  0,  260.00000,  0.00, RGROUP),
   };

   /**
    * Computes the number of implicit hydrogens attached to an atom of type
    * symbol with nsingle single bonds, naromatic aromatic bonds, ndouble
    * double bonds, ntriple trible bonds, and radical and charge state
    * radical and charge, resp.
    */
   public static int implicitHydrogens(String symbol,
                                       int    nsingle,
                                       int    naromatic,
                                       int    ndouble,
                                       int    ntriple,
                                       int    radical,
                                       int    charge)
   {
      int j,h;
      int bond_electrons;

      PTable pt = null;
      for (int i=0; i<ptable.length; i++)
         if (symbol.equals(ptable[i].symbol))
	   {
	      pt = ptable[i];
	      break;
	   }
	   if (pt == null  ||  pt.allowed_valences == null) return (0);

      bond_electrons = nsingle+2*ndouble+3*ntriple;
      if      (radical == Atom.DOUBLET) bond_electrons++;
      else if (radical == Atom.SINGLET) bond_electrons += 2;
      else if (radical == Atom.TRIPLET) bond_electrons += 2;
      switch (naromatic)
      {
         case 0: break;
         case 1: bond_electrons+=2;
                 break;
         case 2: bond_electrons+=3;
                 break;
         case 3: bond_electrons+=4;
                 break;
         default:bond_electrons += naromatic+1;
                 break;
      }

      for (int i=0; i<pt.allowed_valences.length; i++)
      {
         h = pt.allowed_valences[i]-bond_electrons+charge;
         if (0 <= h) return (h);
      }
      return(0);
   }

   private static String HC_table[] =        /* pseudosymbol "G" */
   {"H", "C"};

   private static String non_metal_hetero_elements[] =     /* pseudosymbol "Q" */
   {
      "He",
      "B", "N", "O", "F", "Ne",
      "Si", "P", "S", "Cl", "Ar",
      "As", "Se", "Br", "Kr",
      "Sb", "Te", "I", "Xe",
      "At",                     /* "Rn", This element must be removed */
   };                           /* because of a trick in utils.c */

   private static String metals[] =          /* pseudosymbol "M" */
   {
      "Li", "Be",
      "Na", "Mg", "Al",
      "K", "Ca", "Sc",
        "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
           "Ga",
      "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
           "In", "Sn",
      "Cs", "Ba", "La",
       "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
           "Tl", "Pb", "Bi", "Po",
      "Fr", "Ra", "Ac",
       "Th", "Pa", "U", "Np", "Pu",
  };

   private static String non_metal_small_solution[] =      /* pseudosymbol "Qs" */
   {
      "H",
      "B",  "C", "N", "O", "F",
           "Si", "P", "S", "Cl",
                     "Se", "Br",
                            "I",
   };

   private static String alkali_metals[] =   /* pseudosymbol "alk" */
   {
      "Li", "Na", "K", "Rb", "Cs", "Fr",
  };

   private static String gr2[] =             /* pseudosymbol "gr2" */
   {
      "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
  };

   private static String gr3[] =             /* pseudosymbol "gr3" */
   {
      "B", "Al", "Ga", "In", "Tl",
  };

   private static String gr4[] =             /* pseudosymbol "gr4" */
   {
      "C", "Si", "Ge", "Sn", "Pb",
  };

   private static String ONS_table[] =       /* pseudosymbol "ONS" or "ons" */
   {"O", "N", "S"};

   private static String on2[] =             /* pseudosymbol "on2" */
   {
      "O", "N", "S", "P", "Se", "Te", "Po",
  };

   private static String halogenes[] =       /* pseudosymbol "X" or "hal" */
   {"F", "Cl", "Br", "I", "At"};

   private static String ha2[] =             /* pseudosymbol "ha2" */
   {"Cl", "Br", "I", "At"};

   private static String transition_metals[] =    /* pseudosymbol "trn" */
   {
      "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
      "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
      "La", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
  };

   private static String tra[] =             /* pseudosymbol "tra" */
   {
      "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
      "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
      "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
  };

   private static String trb[] =             /* pseudosymbol "trb" */
   {
      "Cu", "Zn", "Ag", "Cd", "Au", "Hg",
  };

   private static String tm1[] =             /* pseudosymbol "tm1" */
   {
      "Cu", "Ag", "Au",
  };

   private static String tm2[] =             /* pseudosymbol "tm2" */
   {
      "Zn", "Cd", "Hg",
  };

   private static String tm3[] =             /* pseudosymbol "tm3" */
   {
      "Sc", "Y", "La",
  };

   private static String tm4[] =             /* pseudosymbol "tm4" */
   {
      "Ti", "Zr", "Hf",
  };

   private static String tm5[] =             /* pseudosymbol "tm5" */
   {
      "V", "Nb", "Ta",
  };

   private static String tm6[] =             /* pseudosymbol "tm6" */
   {
      "Cr", "Mo", "W",
  };

   private static String tm7[] =             /* pseudosymbol "tm7" */
   {
      "Mn", "Tc", "Re",
  };

   private static String tm8[] =             /* pseudosymbol "tm8" */
   {
      "Fe", "Co", "Ni",
      "Ru", "Rh", "Pd",
      "Os", "Ir", "Pt",
  };

   private static String lanthanoids[] =     /* pseudosymbol "lan" */
   {
       "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
  };

   private static String amino_acids[] =     /* pseudosymbol "Ami" or "ami"*/
   {
      "Ala", "Arg", "Asn", "Asp", "Cys",
      "Gln", "Glu", "Gly", "His", "Ile",
      "Leu", "Lys", "Met", "Phe", "Pro",
      "Ser", "Thr", "Trp", "Tyr", "Val",
   };

   /**
    * Check if sym is equal to one of the array elements of table[].
    */
   private static boolean isInStringTable(String sym, String[] table)
   {
      for (int i=0; i<table.length; i++)
         if (sym.equals(table[i])) return (true);
      return (false);
   }

   /**
    * Returns true if atsym is in the komma delimited list of atom symbols
    * stored in pattern and FALSE otherwise.
    */
   public static boolean atomSymbolMatch(String atsym, String pattern)
   {
      StringTokenizer st = new StringTokenizer(pattern, ",");

      while (st.hasMoreTokens())
      {
         String tok = st.nextToken();
         if (Character.isLowerCase(tok.charAt(0)))
         {
            if ("alk".equals(tok))
            {
               if (isInStringTable(atsym,alkali_metals)) return (true);
            }
            else if ("gr2".equals(tok))
            {
               if (isInStringTable(atsym,gr2)) return (true);
            }
            else if ("gr3".equals(tok))
            {
               if (isInStringTable(atsym,gr3)) return (true);
            }
            else if ("gr4".equals(tok))
            {
               if (isInStringTable(atsym,gr4)) return (true);
            }
            else if ("ons".equals(tok))
            {
               if (isInStringTable(atsym,ONS_table)) return (true);
            }
            else if ("on2".equals(tok))
            {
               if (isInStringTable(atsym,on2)) return (true);
            }
            else if ("hal".equals(tok))
            {
               if (isInStringTable(atsym,halogenes)) return (true);
            }
            else if ("ha2".equals(tok))
            {
               if (isInStringTable(atsym,ha2)) return (true);
            }
            else if ("trn".equals(tok))
            {
               if (isInStringTable(atsym,transition_metals)) return (true);
            }
            else if ("tra".equals(tok))
            {
               if (isInStringTable(atsym,tra)) return (true);
            }
            else if ("trb".equals(tok))
            {
               if (isInStringTable(atsym,trb)) return (true);
            }
            else if ("tm1".equals(tok))
            {
               if (isInStringTable(atsym,tm1)) return (true);
            }
            else if ("tm2".equals(tok))
            {
               if (isInStringTable(atsym,tm2)) return (true);
            }
            else if ("tm3".equals(tok))
            {
               if (isInStringTable(atsym,tm3)) return (true);
            }
            else if ("tm4".equals(tok))
            {
               if (isInStringTable(atsym,tm4)) return (true);
            }
            else if ("tm5".equals(tok))
            {
               if (isInStringTable(atsym,tm5)) return (true);
            }
            else if ("tm6".equals(tok))
            {
               if (isInStringTable(atsym,tm6)) return (true);
            }
            else if ("tm7".equals(tok))
            {
               if (isInStringTable(atsym,tm7)) return (true);
            }
            else if ("tm8".equals(tok))
            {
               if (isInStringTable(atsym,tm8)) return (true);
            }
            else if ("lan".equals(tok))
            {
               if (isInStringTable(atsym,lanthanoids)) return (true);
            }
            else if ("ami".equals(tok))
            {
               if (isInStringTable(atsym,amino_acids)) return (true);
            }
         }
         if (atsym.equals(tok)) return (true);
      }

      if ("A".equals(pattern))
         return (!"H".equals(atsym));
      else if ("Qs".equals(pattern))
         return (isInStringTable(atsym,non_metal_small_solution));
      else if ("G".equals(pattern))
         return (isInStringTable(atsym,HC_table));
      else if ("ONS".equals(pattern))
         return (isInStringTable(atsym,ONS_table));
      else if ("X".equals(pattern))
         return (isInStringTable(atsym,halogenes));
      else if ("Q".equals(pattern))
         return (isInStringTable(atsym,non_metal_hetero_elements));
      else if ("M".equals(pattern))
         return (isInStringTable(atsym,metals));
      else if ("Ami".equals(pattern))
         return (isInStringTable(atsym,amino_acids));
      return (false);
   }
}
