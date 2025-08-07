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

import java.util.Stack;
import java.util.Vector;

import novartis.utilities.Tokenizer;

/**
 * Static functions used to parse a SMILES into a Molecule data structure.
 * Might make it into the Molecule class in the end.
 */
public class SmilesParser
{
   final static int IS_AROMATIC = 1; // flags aromatic atoms in their color field

   /**
    * Helper function to append an integer val to the end of arr[] and return the resulting array.
    */
   private static int[] appendInt(int[] arr, int val)
   {
      int result[] = new int[arr.length+1];
      for (int i=0; i<arr.length; i++) result[i] = arr[i];
      result[result.length-1] = val;
      return (result);
   }

   /**
    * Parse smiles into a Molecule.
    */
   public static Molecule smilesToMolecule(String smiles)
   {
       return smilesToMolecule(smiles, true);
   }

   public static Molecule smilesToMolecule(String smiles, boolean denormalize)
   {
      Vector<Atom> atoms              = new Vector<Atom>();  // temporary storage of Atom objects
      Vector<Bond> bonds              = new Vector<Bond>();  // temporary storage of Bond objects
      int neighbour_numbers[][] = new int[smiles.length()][];  // temporary storage for numbers (index+1) of neighbour atoms
      Stack<Atom>  branch_root_atoms  = new Stack<Atom>();

      Atom open_atoms[]      = new Atom[100];  // only 2 digits for ring connections
      int  open_neighbours[] = new int[100];
      int  open_bond_type[]  = new int[100];

      Tokenizer t = new Tokenizer(smiles);
      String special_atoms[] = {"Cl", "Br", "C", "B", "N", "O", "P", "S", "F", "I",
                                "*",
                                "c", "n", "o", "s", "p"};
      int normal_flags = Tokenizer.IDENT          | Tokenizer.SPECIAL_IDENT |
                         Tokenizer.RING_NUMBER    | Tokenizer.ATTACHMENT    |
                         Tokenizer.BRACKET_STRING | Tokenizer.CHARACTER;

      int bond_type = 0;
      int bond_dir1 = Bond.NORMAL, bond_dir2 = Bond.NORMAL;
      Atom last_atom = null;
      while_loop:
      while (true)
      {
         switch (t.nextToken(normal_flags, special_atoms))
         {
            case Tokenizer.EOF:           // done with SMILES
               break while_loop;
            case Tokenizer.ERROR:         // error in SMILES
               throw (new IllegalArgumentException("Syntax error in SMILES '" +smiles+"'"));
            case Tokenizer.SPECIAL_IDENT:
            case Tokenizer.BRACKET_STRING:
            case Tokenizer.IDENT:
               Atom a;
               if (t.ttype == Tokenizer.SPECIAL_IDENT)         // aromatic atom symbol
               {
                  if (Character.isLowerCase(t.sval.charAt(0)))
                  {
                     StringBuffer symbol = new StringBuffer();
                     symbol.append(t.sval.toUpperCase());
                     if (t.sval.equals("c"))
                        a = new Atom("C", 0, 0, 0, 0);
                     else
                        a = new Atom(0.0, 0.0, 0.0,
                                     symbol.toString(),
                                     null, false,
                                     0,
                                     0, 0, 0, 0);
                     a.color |= IS_AROMATIC;
                  }
                  else if (t.sval.equals("*"))
                  {
                     a = new Atom(0.0, 0.0, 0.0,
                                  "R#", null, false, 0,
                                  0, 0, 0, 0);
                     a.setIntProperty("RGP", 0);
                  }
                  else
                     a = new Atom(t.sval, 0, 0, 0, 0);
               }
               else if (t.ttype == Tokenizer.BRACKET_STRING)   // composed atom symbol
               {
                  a = parseSmilesAtom(t.sval);
               }
               else                                            // normal atom symbol
               {
                  a = new Atom(t.sval, 0, 0, 0, 0);
               }
               atoms.addElement(a);
               neighbour_numbers[atoms.size()-1] = new int[0];
               if (last_atom != null)
               {
                  int at1 = atoms.lastIndexOf(last_atom)+1;
                  int at2 = atoms.size();
                  neighbour_numbers[at1-1] = appendInt(neighbour_numbers[at1-1], at2);
                  neighbour_numbers[at2-1] = appendInt(neighbour_numbers[at2-1], at1);
                  if (bond_type == 0)
                  {
                     if ((((Atom)atoms.elementAt(at1-1)).color & IS_AROMATIC) != 0  &&
                         (((Atom)atoms.elementAt(at2-1)).color & IS_AROMATIC) != 0)
                        bond_type = Bond.AROMATIC;
                     else
                        bond_type = Bond.SINGLE;
                  }
                  Bond b = new Bond(at1, at2, bond_type, Bond.NORMAL, Bond.UNDEFINED, Bond.UNDEFINED);
                  b.color = bond_dir1 + 0x100*bond_dir2;
                  bonds.addElement(b);
               }
               bond_type = 0; bond_dir1 = bond_dir2 = Bond.NORMAL;
               if (a.implicit_H_count == 1)
                  neighbour_numbers[atoms.size()-1] =
                     appendInt(neighbour_numbers[atoms.size()-1], StereoDescription.IMPLICIT_H);
               last_atom = a;
               break;
            case Tokenizer.ATTACHMENT:
               if (last_atom != null)
               {
                  a = new Atom("R#", 0, 0, 0, 0);
                  a.setIntProperty("RGP", t.ival);
                  atoms.addElement(a);
                  neighbour_numbers[atoms.size()-1] = new int[0];
                  int at1 = atoms.lastIndexOf(last_atom)+1;
                  int at2 = atoms.size();
                  neighbour_numbers[at1-1] = appendInt(neighbour_numbers[at1-1], at2);
                  neighbour_numbers[at2-1] = appendInt(neighbour_numbers[at2-1], at1);
                  if (bond_type == 0)  bond_type = Bond.SINGLE;
                  Bond b = new Bond(at1, at2, bond_type, Bond.NORMAL, Bond.UNDEFINED, Bond.UNDEFINED);
                  b.color = bond_dir1 + 0x100*bond_dir2;
                  bonds.addElement(b);
                  bond_type = 0; bond_dir1 = bond_dir2 = Bond.NORMAL;
                  break;
               }
               else
                  throw (new IllegalArgumentException("Illegal attachment in SMILES '" +
                                                      smiles + "'"));
            case Tokenizer.RING_NUMBER:   // opening or closing ring connection
               if (open_atoms[t.ival] == null)  // open connection
               {
                  int at1 = atoms.lastIndexOf(last_atom)+1;
                  neighbour_numbers[at1-1] = appendInt(neighbour_numbers[at1-1], 0);
                  open_atoms[t.ival] = last_atom;
                  open_neighbours[t.ival] = neighbour_numbers[at1-1].length-1;
                  open_bond_type[t.ival] = bond_type;
               }
               else                             // close connection
               {
                  int at1 = atoms.lastIndexOf(open_atoms[t.ival])+1;
                  int at2 = atoms.lastIndexOf(last_atom)+1;
                  neighbour_numbers[at1-1][open_neighbours[t.ival]] = at2;
                  neighbour_numbers[at2-1] = appendInt(neighbour_numbers[at2-1], at1);
                  if (open_bond_type[t.ival] != 0)
                     bond_type = Math.max(bond_type, open_bond_type[t.ival]);
                  if (bond_type == 0)
                  {
                     if ((((Atom)atoms.elementAt(at1-1)).color & IS_AROMATIC) != 0  &&
                         (((Atom)atoms.elementAt(at2-1)).color & IS_AROMATIC) != 0)
                        bond_type = Bond.AROMATIC;
                     else
                        bond_type = Bond.SINGLE;
                  }
                  Bond b = new Bond(at1, at2, bond_type, Bond.NORMAL, Bond.UNDEFINED, Bond.UNDEFINED);
                  b.color = bond_dir1 + 0x100*bond_dir2;
                  bonds.addElement(b);
                  open_atoms[t.ival] = null; open_bond_type[t.ival] = 0;
               }
               bond_type = 0; bond_dir1 = bond_dir2 = Bond.NORMAL;
               break;
            case Tokenizer.CHARACTER:     // bond characters
               switch (t.sval.charAt(0))
               {
                  case '(':   // open branch
                     branch_root_atoms.push(last_atom);
                     bond_type = 0; bond_dir1 = bond_dir2 = Bond.NORMAL;
                     break;
                  case ')':   // close branch
                     last_atom = (Atom)branch_root_atoms.pop();
                     bond_type = 0; bond_dir1 = bond_dir2 = Bond.NORMAL;
                     break;
                  case '-':   // single bond
                     bond_type = Bond.SINGLE;
                     bond_dir1 = bond_dir2; bond_dir2 = Bond.NORMAL;
                     break;
                  case '/':   // single bond UP
                     bond_type = Bond.SINGLE;
                     bond_dir1 = bond_dir2; bond_dir2 = Bond.UP;
                     break;
                  case '\\':  // single bond DOWN
                     bond_type = Bond.SINGLE;
                     bond_dir1 = bond_dir2; bond_dir2 = Bond.DOWN;
                     break;
                  case ':':   // aromatic bond
                     bond_type = Bond.AROMATIC;
                     bond_dir1 = bond_dir2 = Bond.NORMAL;
                     break;
                  case '=':   // double bond
                     bond_type = Bond.DOUBLE;
                     bond_dir1 = bond_dir2 = Bond.NORMAL;
                     break;
                  case '#':   // triple bond
                     bond_type = Bond.TRIPLE;
                     bond_dir1 = bond_dir2 = Bond.NORMAL;
                     break;
                  case '.':   // NOBOND bond
                     bond_type = 0;
                     bond_dir1 = bond_dir2 = Bond.NORMAL;
                     last_atom = null;
                     break;
                  default:
                     throw (new IllegalArgumentException("Illegal character '" +
                                                         t.sval.charAt(0) +
                                                         "' in SMILES '" +
                                                         smiles + "'"));
               }
         }
      }
      // Check if branch stack is empty and bond_type == NOBOND;
      if (!branch_root_atoms.empty())
         throw (new IllegalArgumentException("Branches left open in SMILES '" +smiles+"'"));
      if (bond_type != 0)
         throw (new IllegalArgumentException("Terminal bond character in SMILES '" +smiles+"'"));
      // test if all ring bonds are used up
      for (int i=0; i<100; i++)
         if (open_atoms[i] != null)
            throw (new IllegalArgumentException("Ring bond left open in SMILES '" +smiles+"'"));

      // convert stereo colors/flags into the real Molecule stereo information

      Molecule result = new Molecule();
      result.atoms = new Atom[atoms.size()];
      for (int i=0; i<result.atoms.length; i++)
         result.atoms[i] = (Atom)atoms.elementAt(i);
      result.bonds = new Bond[bonds.size()];
      for (int i=0; i<result.bonds.length; i++)
         result.bonds[i] = (Bond)bonds.elementAt(i);

      decodeSmilesAtomStereo(result, neighbour_numbers);
      decodeSmilesBondStereo(result);

      result.resetColors();

      // Denormalize the molecule, i.e. make aromatic bonds alternating single/double.
      if (denormalize) result.denormalize();

      return (result);
   }

   /**
    * Parses str as a composed SMILES atom and returns a corresponding object.
    */
   static Atom parseSmilesAtom(String str)
   {
      int mass = Atom.NATURAL;
      String symbol = "C";
      int hcount = 0;
      int charge = 0;
      int mapping = Atom.NO_MAPPING;
      String stereo_class = "";
      int    stereo_coset = 0;
      boolean aromatic = false;

      Tokenizer t = new Tokenizer(str);
      String aromatic_atoms[] = {"c", "n", "o", "s", "p"};
      int atom_flags;
      // atom_flags = Tokenizer.IDENT | Tokenizer.SPECIAL_IDENT | Tokenizer.UNSIGNED  | Tokenizer.CHARACTER;

      // expecting mass, atom symbol, or special atom symbol here
      atom_flags = Tokenizer.UNSIGNED | Tokenizer.IDENT | Tokenizer.SPECIAL_IDENT | Tokenizer.CHARACTER;
      t.nextToken(atom_flags, aromatic_atoms);
      if (t.ttype == Tokenizer.UNSIGNED)  // isotopic mass
      {
         mass = t.ival;
         // expecting atom symbol, or special atom symbol here
         atom_flags = Tokenizer.IDENT | Tokenizer.SPECIAL_IDENT | Tokenizer.CHARACTER;
         t.nextToken(atom_flags, aromatic_atoms);
      }

      if (t.ttype == Tokenizer.IDENT          ||
          t.ttype == Tokenizer.SPECIAL_IDENT  ||
          (t.ttype == Tokenizer.CHARACTER  &&  t.sval.equals("*")))
      {
         if (t.ttype == Tokenizer.IDENT  || (t.ttype == Tokenizer.CHARACTER  &&  t.sval.equals("*")))
         {
            symbol = t.sval;
            aromatic = false;
         }
         else
         {
            symbol = t.sval.toUpperCase();
            aromatic = true;
         }
         // expecting stereodesignator ('@'), hcount expression ('H'), charge expression ('+<number>', '-<number>'), or mapping (':<number>') here
         atom_flags = Tokenizer.CHARACTER  | Tokenizer.UNSIGNED;
         t.nextToken(atom_flags, null);
      }
      else
         throw (new IllegalArgumentException("Syntax error in SMILES atom [" +str+"]"));

      if (t.ttype == Tokenizer.CHARACTER  &&  t.sval.equals("@"))          // stereo description
      {
         stereo_class = "TH";
         stereo_coset = 1;
         // expecting stereoclass (e.g. 'TH'), stereocoset ('2' or '@'), or hcount expression ('H') here
         atom_flags = Tokenizer.CHARACTER  | Tokenizer.UNSIGNED;
         t.nextToken(Tokenizer.UCIDENT | Tokenizer.UNSIGNED | Tokenizer.CHARACTER, null);
         if (t.ttype == Tokenizer.CHARACTER  &&   t.sval.equals("@"))      // second tetrahedral coset
         {
            stereo_coset = 2;
            // expecting hcount expression ('H'), charge expression ('+<number>', '-<number>'), or mapping (':<number>') here
            atom_flags = Tokenizer.CHARACTER  | Tokenizer.UNSIGNED;
            t.nextToken(atom_flags, null);
         }
         else if (t.ttype == Tokenizer.UCIDENT  &&  !t.sval.equals("H"))   // class name
         {
            stereo_class = t.sval;
            // expecting hcount expression ('H'), charge expression ('+<number>', '-<number>'), or mapping (':<number>') here
            atom_flags = Tokenizer.CHARACTER;
            t.nextToken(atom_flags, null);
            if (t.ttype == Tokenizer.UNSIGNED)
            {
               stereo_coset = t.ival;
               // expecting hcount expression ('H'), charge expression ('+<number>', '-<number>'), or mapping (':<number>') here
               atom_flags = Tokenizer.CHARACTER;
               t.nextToken(atom_flags, null);
            }
         }
      }

      if (t.ttype == Tokenizer.IDENT   ||
          t.ttype == Tokenizer.UCIDENT ||
          (t.ttype == Tokenizer.CHARACTER  &&  t.sval.equals("H")))
      {
         if (!t.sval.equals("H"))
            throw (new IllegalArgumentException("Syntax error in SMILES atom [" +str+"]"));
         hcount = 1;
         // expecting hcount ('3'), charge expression ('+<number>', '-<number>'), or mapping (':<number>') here
         atom_flags = Tokenizer.CHARACTER  | Tokenizer.UNSIGNED;
         t.nextToken(atom_flags, null);
         if (t.ttype == Tokenizer.UNSIGNED)
         {
            hcount = t.ival;
            t.nextToken(atom_flags, null);
         }
      }

      if (t.ttype == Tokenizer.CHARACTER  &&          // charge description
          (t.sval.equals("+") || t.sval.equals("-")))
      {
         if (t.sval.equals("+"))  charge = 1;
         else                     charge = -1;
         // expecting charge value ('<number>'), or mapping (':<number>') here
         t.nextToken(Tokenizer.UNSIGNED  | Tokenizer.CHARACTER, null);
         if (t.ttype == Tokenizer.UNSIGNED)
         {
            charge *= t.ival;
            // expecting mapping (':<number>') here
            t.nextToken(Tokenizer.CHARACTER, null);
         }
      }

      if (t.ttype == Tokenizer.CHARACTER  &&  t.sval.equals(":")) //mapping
      {
         t.nextToken(Tokenizer.UNSIGNED, null);
         if (t.ttype == Tokenizer.UNSIGNED)
         {
            mapping = t.ival+1;
            t.nextToken(Tokenizer.CHARACTER, null);
         }
         else
            throw (new IllegalArgumentException("Syntax error in SMILES atom [" +str+"]"));
      }

      if (t.ttype != Tokenizer.EOF)
         throw (new IllegalArgumentException("Syntax error in SMILES atom [" +str+"]"));

      Atom result;
      if (symbol.equals("*"))
      {
         result = new Atom(0.0, 0.0, 0.0, "R#",
                           null, false,
                           hcount, mass, 0, Atom.NO_RADICAL, mapping);
         result.setIntProperty("RGP", charge);
      }
      else
         result = new Atom(0.0, 0.0, 0.0, symbol,
                           null, false,
                           hcount, mass, charge, Atom.NO_RADICAL, mapping);

      if (aromatic) result.color |= IS_AROMATIC;
      if (!stereo_class.equals(""))
         result.setIntProperty(stereo_class, stereo_coset);

      return (result);
   }

   /**
    * Convert stereo class atom properties of m to class Molecule internal representation of stereochemistry.
    */
   static void decodeSmilesAtomStereo(Molecule m, int neighbour_lists[][])
   {
      Vector<StereoDescription> result = new Vector<StereoDescription>();

      for (int i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         if (a.getIntProperty("TH", 0) != 0)
         {
            int stereo_coset = a.getIntProperty("TH", 0);
            a.removeIntProperty("TH");
            if (stereo_coset != 1  && stereo_coset != 2)
            {
               System.err.println("Illegal TH coset label '" + stereo_coset + "' ignored");
               continue;
            }

            if (neighbour_lists[i].length != 4)
               System.err.println("tetrahedral stereochemistry for atom with " + neighbour_lists[i].length + " ligands");
            else
            {
               if (stereo_coset == 1)  // listed counter clockwise => swap last two nodes
               {
                  int tmp = neighbour_lists[i][2];
                  neighbour_lists[i][2] = neighbour_lists[i][3];
                  neighbour_lists[i][3] = tmp;
               }
               result.addElement(new StereoDescription(m, a, neighbour_lists[i], StereoDescription.TH));
            }
         }
      }
      if (result.size() > 0)
      {
         m.stereodescs = new StereoDescription[result.size()];
         for (int i=0; i<result.size(); i++)
         m.stereodescs[i] = (StereoDescription)result.elementAt(i);
      }
   }

   /**
    * Convert bond direction colors of m to class Molecule internal representation of stereochemistry.
    */
   static void decodeSmilesBondStereo(Molecule m)
   {
   }

   /**
    * Convert the Bond bond_type to the corresponding string used in SMILES/SMARTS.
    */
   public static String bondTypeToString(int bond_type)
   {
      if (bond_type == 0) return (".");
      if (bond_type == Bond.NOBOND) return (".");

      StringBuffer result = new StringBuffer();
      if ((bond_type&Bond.NOBOND) != 0) result.append('.');
      if ((bond_type&Bond.SINGLE) != 0)
      {
         if (result.length() != 0) result.append(",-");
         else                      result.append('-');
      }
      if ((bond_type&Bond.DOUBLE) != 0)
      {
         if (result.length() != 0) result.append(",=");
         else                      result.append('=');
      }
      if ((bond_type&Bond.TRIPLE) != 0)
      {
         if (result.length() != 0) result.append(",#");
         else                      result.append('#');
      }
      if ((bond_type&Bond.AROMATIC) != 0)
      {
         if (result.length() != 0) result.append(",:");
         else                      result.append(':');
      }

      return (result.toString());
   }

   public static void main(String argv[])
   {
      //Molecule m = smilesToMolecule("c1ccccc1[C@TH1H][NH2]");
      Molecule m = smilesToMolecule("CS(=O)c1c(N)nc(Cl)nc1Cl.CS(=O)c1c(Cl)nc(N)nc1Cl");
   }
}
