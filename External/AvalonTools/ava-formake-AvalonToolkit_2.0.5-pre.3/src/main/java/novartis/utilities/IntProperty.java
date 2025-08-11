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
package novartis.utilities;

/**
 * Represents a integer property list element for an atom (or bond?).
 */
public class IntProperty
{
   public IntProperty next  = null;    // pointer to next element in list
   String             name  = null;    // name of property
   int                value = 0;

   public IntProperty(String name, int value)
   {
      this.name = name;
      this.value = value;
      this.next = null;
   }

   public IntProperty(IntProperty ip)
   {
      this.name = ip.name;
      this.value = ip.value;
      this.next = null;
   }

   public static int findValue(IntProperty prop_list, String name, int default_value)
   {
      while (prop_list != null)
         if (name.equals(prop_list.name)) return (prop_list.value);
         else prop_list = prop_list.next;

      return (default_value);
   }

   /**
    * Sets the property name to value in list.
    * If the named property isn't
    * present, setPropertyInList prepends a new IntProperty object. If the property
    * is on the list, then its value is changed. The function returns the modified
    * property list.
    */
   public static IntProperty setPropertyInList(IntProperty list, String name, int value)
   {
      for (IntProperty ip = list; ip != null; ip=ip.next)
      {
         if (ip.name.equals(name))
         {
            ip.value = value;
            return (list);
         }
      }

      IntProperty ip = new IntProperty(name, value);
      ip.next = list;
      return (ip);
   }

   /**
    * Removes the property called name from list.
    * This is a NOP if the named property isn't present.
    */
   public static IntProperty removePropertyFromList(IntProperty list, String name)
   {
      IntProperty result = null;

      while (list != null)
      {
         if (list.name.equals(name))
            list = list.next;
         else
         {
            IntProperty ip = list.next;
            list.next = result; result = list; list = ip;
         }
      }

      return (result);
   }
}
