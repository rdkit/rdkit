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

import java.util.HashMap;

import novartis.utilities.IntProperty;

/**
 * Represents an atom or general node of a connection table.
 */
public class Atom
{
   double x;
   double y;
   double z;
  

   int    atomic_number;   // contains the standard atomic ordinal number,
                           // if this is a normal atom and some special
                           // ids otherwise.

   String symbol;          // contains standard element symbol from periodic
                           // table if this is a standard element. Special
                           // strings are put here for special node types.
   public String getSymbol() { return symbol; }

   String[] atomList = null; // list of atom symbols used for query atoms
   boolean notLogic = false; // interpretation of list 'false' => 'not in list'
   int isotope;            // isotope of this atom or 0 if natural mix.
   public static final int NATURAL  = 0;
   int charge;             // localized charge of valence bond model.
   /**
    * Get the charge defined by the MDL charge_radical value.
    */
   public static int MDLChargeRadicalToCharge(int charge_radical)
   {
      if      (charge_radical == 1) return (+3);
      else if (charge_radical == 2) return (+2);
      else if (charge_radical == 3) return (+1);
      else if (charge_radical == 5) return (-1);
      else if (charge_radical == 6) return (-2);
      else if (charge_radical == 7) return (-3);
      else                          return (0);
   }

   int radical;            // localized radical state of valence bond model
   // radical state classes. Not identical to MDL's radical types.
   public static final int NO_RADICAL  = 0;
   public static final int SINGLET     = 1;
   public static final int DOUBLET     = 2;
   public static final int TRIPLET     = 4;
   public static final int ANY_RADICAL = 7;
   /**
    * Get the radical state defined by the MDL charge_radical value.
    */
   public static int MDLChargeRadicalToRadical(int charge_radical)
   {
      if (charge_radical == 4) return (DOUBLET);
      else                     return (NO_RADICAL);
   }

   public static final int NO_MAPPING = 0;
   int reaction_mapping = 0;   // reaction map class of this atom if any.
                               // map classes start at 1.

   // list of integer properties of this atom
   IntProperty int_properties = null;
   /**
    * Set an integer property value for this atom.
    */
   public void setIntProperty(String name, int value)
   {
      int_properties =
        IntProperty.setPropertyInList(int_properties, name, value);
   }

   /**
    * Remove an integer property value for this atom.
    */
   public void removeIntProperty(String name)
   {
      int_properties = IntProperty.removePropertyFromList(int_properties, name);
   }

   /**
    * Lookup an integer property value for this atom.
    */
   public int getIntProperty(String name, int default_value)
   {
      return (IntProperty.findValue(int_properties, name, default_value));
   }

   /**
    * Store string keyed properties of atom.
    */
   HashMap<String,String> stringProperties = new HashMap<String,String>();

   /**
    * Setter method for string properties.
    */
   public void setStringProperty(String name, String value)
   {
       stringProperties.put(name, value);
   }

   /**
    * Getter method for string properties.
    */
   public String getStringProperty(String name)
   {
       return (String)stringProperties.get(name);
   }

   // computed fields that are not alway valid.
   final static int DEFAULT_HCOUNT = (-1);
   int implicit_H_count = DEFAULT_HCOUNT;    // number of implict hydrogens attached to this atom.

   public final static int UNDEFINED  = 0; // no topography defined
   public final static int RING       = 1; // atom is defined as a ring atom
   public final static int CHAIN      = 2; // atom is defined as a chain atom
   int topography = UNDEFINED;
  


   // These fields are filled by Molecule.setupNeighbourhood()
   Atom[] neighbour_atoms = null;
   Bond[] neighbour_bonds = null;
   int index = -1;


   // convenience fields for algorithm operating on molecules.
   float value  = 0;
   int   color  = 0;
   int   flags  = 0;
   int   number = 0;

   /**
    * Construct an atom from its common attributes.
    */
   public Atom(double x, double y, double z, String symbol,
               String atomList[], boolean notLogic,
               int implict_H_count,
               int isotope, int charge, int radical, int mapping)
   {
      this.x = x;
      this.y = y;
      this.z = z;
      this.symbol           = symbol;
      this.atomList         = atomList;
      this.notLogic         = notLogic;
      this.implicit_H_count = implict_H_count;
      this.atomic_number    = PTable.SymbolToAtomicNumber(symbol);
      this.isotope          = isotope;
      this.charge           = charge;
      this.radical          = radical;
      this.reaction_mapping = mapping;
   }

   /**
    * Construct an atom from its common attributes w/o coordinates.
    */
   public Atom(String symbol, int isotope, int charge, int radical, int mapping)
   {
      this(0, 0, 0, symbol, null, false,
           DEFAULT_HCOUNT, isotope, charge, radical, mapping);
   }

   /**
    * Deep copy of this object.
    */
   public Object clone()
   {
      // Atom result = new Atom(x, y, z, symbol, implicit_H_count, isotope, charge, radical, reaction_mapping);
      Atom result = new Atom(x, y, z, symbol, atomList, notLogic,
                             implicit_H_count, isotope, charge, radical,
                             reaction_mapping);

      result.value  = value;
      result.color  = color;
      result.flags  = flags;
      result.number = number;

      // clone integer properties
      for (IntProperty ip = int_properties; ip != null; ip = ip.next)
      {
         IntProperty tmp = new IntProperty(ip);
         tmp.next = result.int_properties;
         result.int_properties = tmp;
      }

      // clone String properties
      HashMap<String,String> stringProperties = new HashMap<String,String>();
      for (String key : this.stringProperties.keySet())
      {
          result.stringProperties.put(key, this.stringProperties.get(key));
      }

      return ((Object) result);
   }
}
