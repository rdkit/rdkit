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

import java.util.Vector;

import novartis.utilities.Box;
import novartis.utilities.IntProperty;
import java.awt.Font;

/**
 * This class encapsulates the molecule depiction functions.
 *
 * Depictions are arrays of Strings containing the corresponding commands
 * in a device independent format.
 *
 * The coordinate system is in twips. 1 twip = 1/1440 inch.
 */
public class MoleculeDepicter
{
   // flag for depiction
   public final static int USE_ATOM_COLORS = 0x001;
   public final static int USE_BOND_COLORS = 0x002;
   public final static int USE_COLORS = USE_ATOM_COLORS | USE_BOND_COLORS;

   public final static int USE_CLUSTAL = 0x004;

   // standard bond length is 0.35 inch == 0.889 cm == 504 twips
   public final static int STDBND   = 504;
   // standard font size is 13 points == 260 twips
   public final static int FONTSIZE = 260;

   public final static double CCBND   = 1.54;

   // Java font names of standard fonts.
   // public final static String HELVETICA    = "Helvetica";
   // public final static String TIMES        = "TimesRoman";
   // public final static String COURIER      = "Courier";
   // public final static String DIALOG       = "Dialog";

   // use logical font names
   public final static String HELVETICA    = Font.SANS_SERIF;
   public final static String TIMES        = Font.SERIF;
   public final static String COURIER      = Font.MONOSPACED;
   public final static String DIALOG       = Font.DIALOG;

   /**
    * Select the font for depiction of atom symbols.
    */
   private static void ChooseFont(Vector<String> commands,
                                  double scale,
                                  String fontname)
   {
      commands.addElement(("f " +
                          fontname + " " +
                          (int)(scale*(FONTSIZE+0.5))).toString());
   }

   // depiction bond symbol definitions
   static final int SNG       = 1;
   static final int BND_RIGHT = 2;
   static final int BND_LEFT  = 4;
   static final int DBL_RIGHT = (1|2);
   static final int DBL_LEFT  = (1|4);
   static final int DBL       = DBL_RIGHT;
   static final int TRP       = (1|2|4);
   static final int MANHATTAN = 8192;

   static final int SYM_DBL   = 8;
   static final int WDG       = 16;
   static final int HSH       = 32;
   static final int UNK       = 64;
   static final int DSH       = 128;
   static final int CHN       = 256;
   static final int RNG       = 512;
   static final int BSD       = 1024;
   static final int BSA       = 2048;
   static final int BDA       = 4096;

   static final int CUT       = 1;
   static final int CUT_SEC   = 2;

   /**
    * Creates the commands to draw a bond from (x1,y1) to (x2,y2) and adds it
    * to commands.
    *
    * It uses the font size d as the separation characteristic, the symbol type,
    * and the shortening flags flags1 and flags2 to determine the appearance.
    */
   private static boolean plotBond(Vector<String> commands,
                                   int x1, int y1, int x2, int y2,
                                   int d, int type, int flags1, int flags2)
   {
      int i, len;
      int dbx, dby;       /* multiple bond spacing */
      double bf = 0.45;
      int dwx, dwy;       /* wedge width */
      double wf = 0.5;
      double dhx, dhy;    /* hash line distance */
      int nhash;
      double hf = 0.50;
      int dcx, dcy;       /* line end cutting */
      double cf = 0.65;
      int ddx, ddy;       /* multiple bond cutting */
      double df = 0.45;
      int wedge[][] = {{0, 0}, {0, 0}, {0, 0}};
      int x1a = x1;
      int y1a = y1;
      int x2a = x2;
      int y2a = y2;

      String buf;

      len = (int)Math.sqrt((double)(x1-x2)*(double)(x1-x2) +
                           (double)(y1-y2)*(double)(y1-y2));
      dbx = (int)(bf*((double)d*(x2-x1))/len);
      dby = (int)(bf*((double)d*(y2-y1))/len);
      dwx = (int)(wf*((double)d*(x2-x1))/len);
      dwy = (int)(wf*((double)d*(y2-y1))/len);
      dcx = (int)(cf*((double)d*(x2-x1))/len);
      dcy = (int)(cf*((double)d*(y2-y1))/len);
      ddx = (int)(df*((double)d*(x2-x1))/len);
      ddy = (int)(df*((double)d*(y2-y1))/len);
      // spacing of Manhattan attachments if any

      boolean colorChanged = false;

      if ((type & MANHATTAN)  !=  0)
      {
         if (y1 > y2)   // swap bond ends to make Manhattan bonds align to top
         {
             int h;
             h = x1; x1 = x2; x2 = h;
             h = y1; y1 = y2; y2 = h;
         }
         if (Math.abs(y1-y2) < STDBND*0.5)
         {
             commands.addElement("l " + (int)(x1) + " " + (int)(y1+(d*cf)) + " " +       (int)(x1) + " " + (int)((y1+(d*cf))+d/2));
             // commands.addElement("l " + (int)(x1) + " " + (int)((y1+(d*cf))+d/2) + " " + (int)(x2) + " " + (int)((y2-(d*cf))-d/2));
             commands.addElement("l " + (int)(x1) + " " + (int)((y1+(d*cf))+d/2) + " " + (int)(x2) + " " + (int)((y2+(d*cf))+d/2));
             // commands.addElement("l " + (int)(x2) + " " + (int)((y2-(d*cf))-d/2) + " " + (int)(x2) + " " + (int)(y2-(d*cf)));
             commands.addElement("l " + (int)(x2) + " " + (int)(y2+(d*cf)) + " " +       (int)(x2) + " " + (int)((y2+(d*cf))+d/2));
         }
         else
         {
             commands.addElement("l " + (int)(x1) + " " + (int)(y1+(d*cf)) + " " +       (int)(x1) + " " + (int)((y1+(d*cf))+d/2));
             commands.addElement("l " + (int)(x1) + " " + (int)((y1+(d*cf))+d/2) + " " + (int)(x2) + " " + (int)((y2-(d*cf))-d/2));
             commands.addElement("l " + (int)(x2) + " " + (int)((y2-(d*cf))-d/2) + " " + (int)(x2) + " " + (int)(y2-(d*cf)));
         }
         return colorChanged;
      }

      if ((flags1 & CUT)  !=  0) { x1 += dcx; y1 += dcy; }
      if ((flags2 & CUT)  !=  0) { x2 -= dcx; y2 -= dcy; }

      if ((type & SNG)  !=  0)
         commands.addElement("l " + x1 + " " + y1 + " " + x2 + " " + y2);

      if ((type & BND_RIGHT)  !=  0)
      {
         if ((flags1 & CUT_SEC)  !=  0)
            buf = "l " + (x1-dby+ddx) + " " + (y1+dbx+ddy) + " ";
         else
            buf = "l " + (x1-dby) + " " + (y1+dbx) + " ";
         if ((flags2 & CUT_SEC)  !=  0)
            buf = buf + (x2-dby-ddx) + " " + (y2+dbx-ddy);
         else
            buf = buf + (x2-dby) + " " + (y2+dbx);
         commands.addElement(buf);
      }

      if ((type & BND_LEFT)  !=  0)
      {
         if ((flags1 & CUT_SEC)  !=  0)
            buf = "l " + (x1+dby+ddx) + " " + (y1-dbx+ddy) + " ";
         else
            buf = "l " + (x1+dby) + " " + (y1-dbx) + " ";
         if ((flags2 & CUT_SEC)  !=  0)
            buf = buf + (x2+dby-ddx) + " " + (y2-dbx-ddy);
         else
            buf = buf + (x2+dby) + " " + (y2-dbx);
         commands.addElement(buf);
      }

      if ((type & SYM_DBL)  !=  0)
      {
         buf = "l " + (x1-dby/2) + " " + (y1+dbx/2) + " " +
                      (x2-dby/2) + " " + (y2+dbx/2);
         commands.addElement(buf);
         buf = "l " + (x1+dby/2) + " " + (y1-dbx/2) + " " +
                      (x2+dby/2) + " " + (y2-dbx/2);
         commands.addElement(buf);
      }

      if ((type & WDG)  !=  0)
      {
         wedge[0][0] = x1; wedge[0][1] = y1;
         wedge[1][0] = x2-dwy/2; wedge[1][1] = y2+dwx/2;
         wedge[2][0] = x2+dwy/2; wedge[2][1] = y2-dwx/2;
         buf = "pf 3 " +
               wedge[0][0] + " " + wedge[0][1] + " " +
               wedge[1][0] + " " + wedge[1][1] + " " +
               wedge[2][0] + " " + wedge[2][1];
         commands.addElement(buf);
      }

      if ((type & HSH)  !=  0)
      {
         nhash = (int)(len/(hf*d)+0.5);
         if (nhash < 2) nhash = 2;
         dhx = wf*((double)d*(x2-x1))/len;
         dhy = wf*((double)d*(y2-y1))/len;
         for (i=0; i<=nhash; i++)
         {
            buf = "l " + (int)(x1+i*(x2-x1-dhy/2)/nhash) + " " +
                         (int)(y1+i*(y2-y1+dhx/2)/nhash) + " " +
                         (int)(x1+i*(x2-x1+dhy/2)/nhash) + " " +
                         (int)(y1+i*(y2-y1-dhx/2)/nhash);
            commands.addElement(buf);
         }
      }

      if ((type & DSH)  !=  0)
      {
         nhash = (int)(len/(hf*d)+0.5);
         if (nhash < 4) nhash = 4;
         if (nhash%2 == 1) nhash++;
         if ((type&SNG) != 0)
         {
             x1 += -dby+ddx; y1 -= -dbx-ddy;
             x2 += -dby-ddx; y2 -= -dbx+ddy;
         }
         for (i=0; i<nhash; i+=2)
         {
            buf = "l " + (int)(x1+i*(x2-x1)/nhash) + " " +
                         (int)(y1+i*(y2-y1)/nhash) + " " +
                         (int)(x1+(i+1)*(x2-x1)/nhash) + " " +
                         (int)(y1+(i+1)*(y2-y1)/nhash);
            commands.addElement(buf);
         }
      }

      if ((type & CHN) != 0)
      {
          colorChanged = true;
          commands.addElement("c " + colortable[1]);    // make label red
          commands.addElement("t 0 " + ((x1a+x2a)/2) + " " +
                                       ((y1a+y2a)/2) + " ch");
      }

      if ((type & RNG) != 0)
      {
          colorChanged = true;
          commands.addElement("c " + colortable[1]);    // make label red
          commands.addElement("t 0 " + ((x1a+x2a)/2) + " " +
                                       ((y1a+y2a)/2) + " rn");
      }

      if ((type & BSD) != 0)
      {
          colorChanged = true;
          commands.addElement("c " + colortable[1]);    // make label red
          commands.addElement("t 0 " + ((x1a+x2a)/2) + " " +
                                       ((y1a+y2a)/2) + " s/d");
      }

      if ((type & BSA) != 0)
      {
          colorChanged = true;
          commands.addElement("c " + colortable[1]);    // make label red
          commands.addElement("t 0 " + ((x1a+x2a)/2) + " " +
                                       ((y1a+y2a)/2) + " s/a");
      }

      if ((type & BDA) != 0)
      {
          colorChanged = true;
          commands.addElement("c " + colortable[1]);    // make label red
          commands.addElement("t 0 " + ((x1a+x2a)/2) + " " +
                                       ((y1a+y2a)/2) + " d/a");
      }
      return colorChanged;
   }

   /*
    * Computes the drawing window of the molecule m in world coordinates.
    * Returns the second decile bond length.
    *
    * The median bond length tends to be too large at times. 
    */
   public static double getWorldWindow(Molecule m, Box box)
   {
      Atom a;
      Bond b;
      int i;
      double d1, d2, len;

      box.xmin = box.ymin = 1e10; box.xmax = box.ymax = -1e10;      /* World Window */
      for (i=0; m.atoms != null  &&  i<m.atoms.length; i++)
      {
         a = m.atoms[i];
         if (box.xmin > a.x) box.xmin = a.x;
         if (box.xmax < a.x) box.xmax = a.x;
         if (box.ymin > a.y) box.ymin = a.y;
         if (box.ymax < a.y) box.ymax = a.y;
      }
      if (box.xmin == box.xmax) {box.xmin -= 0.5; box.xmax += 0.5;}
      if (box.ymin == box.ymax) {box.ymin -= 0.5; box.ymax += 0.5;}

      len = 0.0;
      double[] lengths = new double[1];
      lengths[0] = CCBND/2.0;   // make sure there is at least one length record (of half a CC bond)
      if (m.bonds != null  &&  m.bonds.length > 0)
          lengths = new double[m.bonds.length];
      int nshortcut = 0;
      for (i=0; m.bonds != null  &&  i<m.bonds.length; i++)
      {
         b = m.bonds[i];
         Atom a1 = m.atoms[b.atoms[0]-1];
         Atom a2 = m.atoms[b.atoms[1]-1];
         boolean isShortcutBond = false;
         String s1 = a1.symbol;
         String s2 = a2.symbol;
         if (a1.getStringProperty("A") != null)
         {
             isShortcutBond = true;
             s1 = a1.getStringProperty("A");
         }
         if (a2.getStringProperty("A") != null)
         {
             isShortcutBond = true;
             s2 = a2.getStringProperty("A");
         }
         d1 = a1.x - a2.x; d2 = a1.y - a2.y;
         lengths[i] = d1*d1 + d2*d2;
         if (isShortcutBond)    // shotcut bonds are elongated for longer labels
         {
             nshortcut++;
             lengths[i] *= 4.0/(3.0+0.5*(s1.length()+s2.length()));
             lengths[i] *= 4.0/(3.0+0.5*(s1.length()+s2.length()));
         }
      }
      for (i=1; i<lengths.length; i++)
          for (int j=i-1; j>=0; j--)
              if (lengths[j+1] < lengths[j])
              {
                  double tmp = lengths[j+1]; lengths[j+1] = lengths[j]; lengths[j] = tmp;
              }
              else
                  break;
      len = Math.sqrt(lengths[2*lengths.length/10]);

      // increase right bound if needed to provide for atom labels
      for (i=0; m.atoms != null  &&  i<m.atoms.length; i++)
      {
         a = m.atoms[i];
         String s = a.symbol;
         if (a.getStringProperty("A") != null)
         {
             s = a.getStringProperty("A");
         }
         else if (a.atomList != null)
         {
             StringBuffer buffer = new StringBuffer();
             for (int j=0; j<a.atomList.length; j++)
             {
                 if (j==0)
                     buffer.append("[");
                 else
                 {
                     if (a.notLogic)
                        buffer.append(";");
                     else
                        buffer.append(",");
                 }
                 if (a.notLogic) buffer.append("!");
                 buffer.append(a.atomList[j]);
             }
             buffer.append("]");
             s = buffer.toString();
         }
         if (box.xmax < a.x+s.length()*len/4.0)
         {
             box.xmax = a.x+s.length()*len/4.0;
         }
      }

      // add some space to left and right of box
      box.xmin -= len/2.0; box.xmax += len/2.0;
      return (len);
   }

   // instance variables
   double xoffset;
   // instance variables
   double yoffset;
   double scale = 1.0 ;
   double fontscale = 1.0 ;

   /**
    * Set coordinate transformation.
    */
   private void setTransformation(int metalen,
                                  Box b,
                                  double len)
   {
      scale = metalen/len;
      fontscale = 1.0;
      xoffset = (b.xmin+b.xmax)/2;
      yoffset = (b.ymin+b.ymax)/2;
   }

   private void fixScale(double xext_best, double yext_best,
                         double xext_new,  double yext_new)
   /*
    * Sets the global scale variable such that the picture will now
    * fit into the box defined by xext_new and yext_new.
    */
   {
      if (xext_best < xext_new  &&  yext_best < yext_new)
      {
         return;     /* best size would already fix into new box => NOOP */
      }

      if (xext_new/xext_best <= yext_new/yext_best)
      {
         /* xscale defines scale */
         scale *= xext_new/xext_best;
         fontscale *= xext_new/xext_best;
      }
      else
      {
         /* yscale defines scale */
         scale *= yext_new/yext_best;
         fontscale *= yext_new/yext_best;
      }
      // if (fontscale < 0.3) fontscale = 0.3;
      // if (fontscale < 0.5) fontscale = 0.5;
      // if (fontscale < 0.8) fontscale = 0.8;
   }

   /**
    * Convert world coordinate X to meta coordinate system.
    */
   int worldToMetaX(double x)
   {
      return ((int)((x-xoffset)*scale));
   }

   /**
    * Convert world coordinate Y to meta coordinate system.
    */
   int worldToMetaY(double y)
   {
      return (-(int)((y-yoffset)*scale));
   }

   static void SetDrawFlags(Molecule m)
   /*
    * This procedure looks for each atom and bond how they should be drawn.
    * The information is stored in the color fields of atoms and bonds.
    */
   {
      Atom a, a1, a2;
      Bond b;
      int at1, at2;

      m.resetColors();
      for (int i=0; i<m.bonds.length; i++)
      {
         b = m.bonds[i];
         b.color = 0;
         if (b.topography == Atom.CHAIN) b.color |= CHN;
         if (b.topography == Atom.RING)  b.color |= RNG;
      }

      m.setupNeighbourhood();
      m.perceiveRingBonds();

      m.makeRingsClockwise();

      for (int i=0; i<m.atoms.length; i++)
      {
         a = m.atoms[i];
         if (!a.symbol.equals("C")                      ||
             a.charge  != 0                             ||
             a.radical != Atom.NO_RADICAL               ||
             (a.implicit_H_count != Atom.DEFAULT_HCOUNT && a.reaction_mapping == 0) ||   // not really correct, but used for first version
             a.getIntProperty("SUB", 0) != 0            ||
             a.isotope != 0)
             a.color |= CUT;
         else if (m.atoms[i].neighbour_atoms.length != 1)
             a.color |= CUT_SEC;
      }

      for (int i=0; i<m.bonds.length; i++)
      {
         b = m.bonds[i];
         switch (b.bond_type)
         {
            case Bond.SINGLE: if (b.stereo_symbol == Bond.UP)
                                 b.color |= WDG;
                              else if (b.stereo_symbol == Bond.DOWN)
                                 b.color |= HSH;
                              else
                                 b.color |= SNG;
                              break;
            case Bond.DOUBLE: at1 = b.atoms[0]; at2 = b.atoms[1];
                              a1 = m.atoms[at1-1];
                              a2 = m.atoms[at2-1];
                              if (b.topography != Bond.RING)
                              {     /* dbl bond should have common direction */
                                 if ((a2.x-a1.x)*2 + (a2.y-a1.y)*3  < 0)
                                 {
                                    b.atoms[0] = at2; b.atoms[1] = at1;
                                 }
                              }
                              if (a1.neighbour_atoms.length < 2   ||
                                  a2.neighbour_atoms.length < 2   ||
                                  (a1.neighbour_atoms.length == 2 && a2.neighbour_atoms.length == 2 &&
                                   (!a1.symbol.equals("C")  || !a2.symbol.equals("C"))  &&
                                   b.topography != Bond.RING)     ||
                                  (a1.neighbour_atoms.length == 3 && a2.neighbour_atoms.length == 3 &&
                                   // a1.topography == Atom.RING     && a2.topography == Atom.RING  &&
                                   b.topography != Bond.RING))
                                 b.color |= SYM_DBL;
                              else
                                 b.color |= DBL;
                              break;
            case Bond.TRIPLE: b.color |= TRP;
                              at1 = b.atoms[0]; at2 = b.atoms[1];
                              a1 = m.atoms[at1-1];
                              a2 = m.atoms[at2-1];
                              a1.color &= ~CUT_SEC; a2.color &= ~CUT_SEC;
                              break;
            case Bond.SINGLE|Bond.DOUBLE:
                              b.color |= BSD;
                              b.color |= DSH;
                              b.color |= SNG;
                              break;
            case Bond.SINGLE|Bond.AROMATIC:
                              b.color |= BSA;
                              b.color |= DSH;
                              b.color |= SNG;
                              break;
            case Bond.DOUBLE|Bond.AROMATIC:
                              b.color |= BDA;
                              b.color |= DSH;
                              b.color |= SNG;
                              break;
            case Bond.AROMATIC:
                              b.color |= DSH;
                              b.color |= SNG;
                              break;
            default: b.color |= DSH;
                     break;
         }
      }
   }

   String AttributesToString(int charge,
                             int radical,
                             int mass_difference)
   {
      StringBuffer result = new StringBuffer();

      if (mass_difference != 0) result.append("*");

      if (Math.abs(charge) > 1) result.append(""+Math.abs(charge));

      if (charge > 0)
         result.append("+");
      else if (charge < 0)
         result.append("-");

      if      (radical == Atom.SINGLET)
         result.append("|");
      else if (radical == Atom.DOUBLET)
         result.append(".");
      else if (radical == Atom.TRIPLET)
         result.append(":");

      return (result.toString());
   }

   /**
    * Lookup table for standard colors.
    */
   static String colortable[] =
      {
         "0x000000",        // black
         "0xFF0000",        // red
         "0x0000FF",        // blue
         "0x007F00",        // dark green
         "0x7F3F00",        // brown
         "0x3F007F",        // dark violet
         "0x003F3F",        // dark cyan
         "0x7F7F00",        // dark yellow
         "0x7F007F",        // dark pink
         "0x999999",        // grey
      };

   static final int GREY_INDEX = 9;

   /**
    * Lookup table for shortcut colors.
    */
   static String clustal_colors[] =
      {
         "0xf09048",        // ORANGE
         "0xc0c000",        // YELLOW
         "0x80a0f0",        // BLUE
         "0xf01505",        // RED
         "0x15c015",        // GREEN
         "0xf08080",        // PINK
         "0xc048c0",        // MAGENTA
         "0x15a4a4",        // CYAN
      };

   /**
    * Lookup table for shortcut colors.
    */
   static String rasmol_colors[] =
      {
         "0xe60a0a",        // bright red
         "0x145aff",        // blue
         "0x3232aa",        // mid blue
         // "0xebebeb",        // light grey
         // "0xc8c8c8",        // dark grey
         // need to make them darker for white background
         "0xc8c8c8",        // light grey
         "0xa4a4a4",        // dark grey
         "0x8282d2",        // pale blue
         // "0xe6e600",        // yellow (using Clustal color for better readability)
         "0xc0c000",        // YELLOW
         "0xfa9600",        // orange
         "0x00dcdc",        // cyan
         "0x0f820f",        // green
         "0xb45ab4",        // pink
         "0xdc96dc",        // flesh
      };

    /**
     * Find the best direction into which to render the atom symbol.
     *
     * It uses the neighbour coordinates to determine where not bonds
     * would be in the way.
     */
    private String freeDirection(Atom a, Atom[] neighbours)
    {
        // First, we look for a completely free direction,
        // e.g. for terminal bonds.
        int j;
        if (a.atomList != null) return "E";
        // regular shortcuts are always left-to-right, but repeat shortcuts can be otherwise
        if (a.getStringProperty("A") != null  &&  !a.getStringProperty("A").matches(".*[0-9].*")) return "E";
        if (a.getStringProperty("A") != null  &&  a.getStringProperty("A").matches(".*[A-Z][A-Z].*")) return "E";
        if (a.getIntProperty("SUB", 0) != 0) return "E";
        for (j=0; j<neighbours.length; j++)
           if (neighbours[j].x-a.x > 0) break;
        if (j >= neighbours.length) return "E";
        for (j=0; j<neighbours.length; j++)
           if (neighbours[j].x-a.x < 0) break;
        if (j >= neighbours.length) return "W";
        for (j=0; j<neighbours.length; j++)
           if (neighbours[j].y-a.y < 0) break;
        if (j >= neighbours.length) return "S";
        for (j=0; j<neighbours.length; j++)
           if (neighbours[j].y-a.y > 0) break;
        if (j >= neighbours.length) return "N";

        // Now, we look for smaller spaces
        for (j=0; j<neighbours.length; j++)
        {
           if (neighbours[j].x-a.x > 0  &&
               Math.abs(neighbours[j].x-a.x) > Math.abs(neighbours[j].y-a.y))
              break;
        }
        if (j >= neighbours.length) return "e";
        for (j=0; j<neighbours.length; j++)
        {
           if (neighbours[j].x-a.x < 0  &&
               Math.abs(neighbours[j].x-a.x) > Math.abs(neighbours[j].y-a.y))
              break;
        }
        if (j >= neighbours.length) return "w";
        for (j=0; j<neighbours.length; j++)
        {
           if (neighbours[j].y-a.y < 0  &&
               Math.abs(neighbours[j].x-a.x) < Math.abs(neighbours[j].y-a.y))
              break;
        }
        if (j >= neighbours.length) return "s";
        for (j=0; j<neighbours.length; j++)
        {
           if (neighbours[j].y-a.y > 0  &&
               Math.abs(neighbours[j].x-a.x) < Math.abs(neighbours[j].y-a.y))
              break;
        }
        if (j >= neighbours.length) return "n";
        return "e";
    }

   /**
    * Chemical structure depiction workhorse. Takes the Molecule m
    * and returns an array of device independent drawing commands.
    *
    * The flags parameter currently tells if color is to be used.
    *
    * The array colors[] is used to color the different parts of the molecule.
    * colors[i] == 0 means draw in black. 1 <= colors[i] <= 9 means use one of
    * the nine most dark colors to get reasonable readablility. If colors ==
    * NULL, the standard atom label coloring with black bonds is used.
    *
    * xext and yext are used for the desired x- and y- dimensions of the resulting
    * image in twips.
    *
    * The array labels[] contains alternate atom labels to be used instead
    * of the ones defined in the symbol fields.
    * Labels which are "" strings are not considered. Use " " to get
    * a blank label.
    */
   public String[] computeDepiction(Molecule m,
                                    int      xext,
                                    int      yext,
                                    int      flags,
                                    int      colors[],
                                    String   labels[])
   {
      Bond b;

      Box box = new Box(0, 0, 0, 0);
      double len;
      boolean color_changed=true;
      int xExt, yExt;
      int x1, y1, x2, y2;
      int flags1, flags2;
      String atom_symbol = "";
      String symbolp;
      int charge;
      int[] H_count;

      int current_color;

      Vector<String> commands;

      if (m == null) return (new String[0]);

      len = getWorldWindow(m, box);

      setTransformation(STDBND, box, len);

      xExt = worldToMetaX(box.xmax)-worldToMetaX(box.xmin);
      yExt = worldToMetaY(box.ymin)-worldToMetaY(box.ymax);
      xExt += 1.5*STDBND; yExt += 1.5*STDBND;

      if ((xext) > 0  &&  (yext) > 0)
      {
         fixScale((double)xExt,   (double)yExt,
                  (double)(xext), (double)(yext));
         xExt = xext; yExt = yext;
      }

      commands = new Vector<String>();

      commands.addElement("e " + xExt + " " + yExt);
      commands.addElement("o " + -xExt/2 + " " + -yExt/2);

      current_color = 0;
      commands.addElement("c " + colortable[current_color % colortable.length]);

      SetDrawFlags(m);

      /* Draw bonds */
      for (int i=0; i<m.bonds.length; i++)
      {
         b = m.bonds[i];
         double blen = 0;
         blen += (m.atoms[b.atoms[0]-1].x-m.atoms[b.atoms[1]-1].x)*(m.atoms[b.atoms[0]-1].x-m.atoms[b.atoms[1]-1].x);
         blen += (m.atoms[b.atoms[0]-1].y-m.atoms[b.atoms[1]-1].y)*(m.atoms[b.atoms[0]-1].y-m.atoms[b.atoms[1]-1].y);
         blen = Math.sqrt(blen);
         x1 = worldToMetaX(m.atoms[b.atoms[0]-1].x);
         y1 = worldToMetaY(m.atoms[b.atoms[0]-1].y);
         x2 = worldToMetaX(m.atoms[b.atoms[1]-1].x);
         y2 = worldToMetaY(m.atoms[b.atoms[1]-1].y);
         flags1 = m.atoms[b.atoms[0]-1].color;
         flags2 = m.atoms[b.atoms[1]-1].color;
         // don't draw horizontal bonds connecting shortcuts except for ring closure bonds spanning more than 1.9 average bond distances
         if (Math.abs(y1-y2) < 0.001*Math.abs(x1-x2)  &&  STDBND*2.9 > Math.abs(x1-x2)  &&  b.bond_type == Bond.SINGLE)
         { 
             boolean isShortcutLink = false;
             Atom a = m.atoms[b.atoms[0]-1];
             if ((a.symbol.equals("R") || a.symbol.equals("A"))  &&  null != a.getStringProperty("A")) isShortcutLink = true;
             a = m.atoms[b.atoms[1]-1];
             if ((a.symbol.equals("R") || a.symbol.equals("A"))  &&  null != a.getStringProperty("A")) isShortcutLink = true;
             // clear bond only if both ends are shortcuts
             a = m.atoms[b.atoms[0]-1];
             if (null == a.getStringProperty("A")) isShortcutLink = false;
             if (!(a.symbol.equals("R") || a.symbol.equals("A")))
             {
                 if (a.symbol.length() > 1) isShortcutLink = false;
             }
             a = m.atoms[b.atoms[1]-1];
             if (null == a.getStringProperty("A")) isShortcutLink = false;
             if (!(a.symbol.equals("R") || a.symbol.equals("A")))
             {
                 if (a.symbol.length() > 1) isShortcutLink = false;
             }
             if (isShortcutLink) continue;
         }
         if (colors != null)    /* use coded colors */
         {
            int new_color = Math.min(colors[b.atoms[0]-1],
                                     colors[b.atoms[1]-1]);
            if (new_color < colortable.length)
            {   /* use this color */
               if (current_color != new_color)
               {
                  current_color = new_color;
                  commands.addElement("c " + colortable[current_color % colortable.length]);
               }
            }
            else
            {   /* use black */
               if (current_color != 0)
               {
                  current_color = 0;
                  commands.addElement("c " +
                  colortable[current_color % colortable.length]);
               }
            }
         }
         else if (blen > len*2)     // draw long bonds linking shortcuts in grey and as Manhattan bonds
         {
             boolean isShortcutLink = true;
             Atom a = m.atoms[b.atoms[0]-1];
             if (!(a.symbol.equals("R") || a.symbol.equals("A"))  ||  null == a.getStringProperty("A")) isShortcutLink = false;
             a = m.atoms[b.atoms[1]-1];
             if (!(a.symbol.equals("R") || a.symbol.equals("A"))  ||  null == a.getStringProperty("A")) isShortcutLink = false;
             if (isShortcutLink)        // both ends are shortcut symbols
             {
                 commands.addElement("c " + colortable[GREY_INDEX]);
                 current_color = GREY_INDEX;
                 b.color |= MANHATTAN;
             }
             else
             {
                 commands.addElement("c " + colortable[0]);
                 current_color = 0;
             }
         }
         else if ((flags & USE_BOND_COLORS) != 0)
         {
            int new_color = 0;
            if (b.reaction_mark > 0) new_color = 1;
            if (current_color != new_color)
            {
               current_color = new_color;
               commands.addElement("c " + colortable[new_color]);
            }
         }
         else
         {
            if (current_color != 0)
            {
               current_color = 0;
               commands.addElement("c " + colortable[current_color % colortable.length]);
            }
         }
         if (plotBond(commands,
                  x1, y1, x2, y2,
                  (int)(fontscale*FONTSIZE), b.color, flags1, flags2))
                  commands.addElement("c " + colortable[current_color]);
      }

      current_color = 0;
      commands.addElement("c " + colortable[current_color]);


      /* draw atom symbol strings */
      ChooseFont(commands, fontscale, HELVETICA);

      H_count = m.ComputeImplicitH();

      current_color = 0;
      commands.addElement("c " + colortable[0]);
      color_changed = false;
      String direction = "";
      for (int i=0; i<m.atoms.length; i++)
      {
         Atom a = m.atoms[i];
         // mark atom coordinates
         commands.addElement("ac " + (i+1) + " " + worldToMetaX(a.x) + " " + worldToMetaY(a.y));
         // print symbol if needed
         if (!a.symbol.equals("C")                                                 ||
             a.charge  != 0                                                        ||
             a.radical != Atom.NO_RADICAL                                          ||
             a.isotope != Atom.NATURAL                                             ||
             m.atoms.length == 1                                                   ||
             (labels    != null  &&  labels[i] != null  &&  !labels[i].equals("")) ||
              a.getIntProperty("SUB", 0) != 0                                      ||
              a.atomList != null)
         {
            if (a.charge  == 0               &&
                a.radical == Atom.NO_RADICAL &&
                a.isotope == Atom.NATURAL    &&
                a.atomList == null)
            {
               if      (a.symbol.equals("Cl") && H_count[i+1] == 1)
                  atom_symbol = "HCl";
               else if (a.symbol.equals("Br") && H_count[i+1] == 1)
                  atom_symbol = "HBr";
               else if (a.symbol.equals("I")  && H_count[i+1] == 1)
                  atom_symbol = "HI";
               else if (a.symbol.equals("F")  && H_count[i+1] == 1)
                  atom_symbol = "HF";
               else if (a.symbol.equals("O")  && H_count[i+1] == 2)
                  atom_symbol = "H2O";
               else if (a.symbol.equals("S")  && H_count[i+1] == 2)
                  atom_symbol = "H2S";
               else if ((a.symbol.equals("R") || a.symbol.equals("A"))  && null != a.getStringProperty("A"))
               {
                  atom_symbol = a.getStringProperty("A");
               }
               else if (H_count[i+1] == 0  ||  a.getIntProperty("SUB", 0) != 0)
                  atom_symbol = a.symbol;
               else if (H_count[i+1] == 1)
                  atom_symbol =  a.symbol + "H";
               else if (H_count[i+1] > 1)
                  atom_symbol = a.symbol + "H" + H_count[i+1];
            }
            else if (a.atomList != null)
            {
                StringBuffer buffer = new StringBuffer();
                for (int j=0; j<a.atomList.length; j++)
                {
                    if (j==0)
                        buffer.append("[");
                    else
                    {
                        if (a.notLogic)
                           buffer.append(";");
                        else
                           buffer.append(",");
                    }
                    if (a.notLogic) buffer.append("!");
                    buffer.append(a.atomList[j]);
                }
                buffer.append("]");
                atom_symbol = buffer.toString();
            }
            else
            {
               if (H_count[i+1] == 0)
                atom_symbol = a.symbol;
               else if (H_count[i+1] == 1)
                atom_symbol = a.symbol + "H";
               else if (H_count[i+1] > 1)
                atom_symbol = a.symbol + "H" + H_count[i+1];
            }
            if (a.getIntProperty("SUB", 0) != 0)
            {
                int property = a.getIntProperty("SUB", 0);
                if (property == -2)
                    atom_symbol = atom_symbol + "(s*)";
                else if (property == -1)
                    atom_symbol = atom_symbol + "(s0)";
                else if (property > 0)
                    atom_symbol = atom_symbol + "(s"+property+")";
            }

            Atom[] neighbours = m.atoms[i].neighbour_atoms;
            // Find a free direction for rendering the symbol
            direction = freeDirection(a, neighbours);

            charge = a.charge;
            if (a.atomic_number == PTable.RGROUP)
            {
               if      (charge == -1)
               {
                atom_symbol = "From";
               }
               else if (charge == -2)
               {
                atom_symbol = "To";
               }
               else if (charge < 0)
               {
                atom_symbol = "R?";
               }
               else
               {
                  charge = IntProperty.findValue(a.int_properties, "RGP", 0);
                  atom_symbol = "R" + charge;
                  charge = 0;
               }
            }
            if ((flags & USE_ATOM_COLORS) != 0  &&  (colors == null))
            {
               color_changed = true;
               if (a.symbol.equals("N"))          /* blue */
               {
                commands.addElement("c " + colortable[2]);
               }
               else if (a.symbol.equals("O"))     /* red */
               {
                commands.addElement("c " + colortable[1]);
               }
               else if (a.symbol.equals("S"))     /* brown */
               {
                commands.addElement("c " + colortable[4]);
               }
               else if (a.symbol.equals("P"))     /* dark pink */
               {
                  commands.addElement("c " + colortable[8]);
               }
               else if (a.symbol.equals("I"))     /* violet */
               {
                  commands.addElement("c " + colortable[5]);
               }
               else if (a.symbol.equals("F"))     /* green */
               {
                  commands.addElement("c " + colortable[3]);
               }
               else if (a.symbol.equals("Cl"))    /* green */
               {
                  commands.addElement("c " + colortable[3]);
               }
               else if (a.symbol.equals("Br"))    /* green */
               {
                  commands.addElement("c " + colortable[3]);
               }
               else if (a.symbol.equals("R")  &&  0 == (flags & USE_CLUSTAL))    /* shortcuts colored the RASMOL way */
               {
                   String shortcut = a.getStringProperty("A");
                   if (shortcut == null) 
                      commands.addElement("c " + colortable[0]);
                   else if (shortcut.equals("E")  ||  shortcut.matches(".*\\b[D]?(Me)?Glu")  ||
                            shortcut.equals("D")  ||  shortcut.matches(".*\\b[D]?(Me)?Asp"))
                      commands.addElement("c " + rasmol_colors[0]);    // bright red
                   else if (shortcut.equals("R")  ||  shortcut.matches(".*\\b[D]?(Me)?Arg")  ||
                            shortcut.equals("K")  ||  shortcut.matches(".*\\b[D]?(Me)?Lys"))
                      commands.addElement("c " + rasmol_colors[1]);    // blue
                   else if (shortcut.equals("F")  ||  shortcut.matches(".*\\b[D]?(Me)?Phe")  ||
                                                      shortcut.matches(".*\\b[D]?(Me)?aNal") ||
                                                      shortcut.matches(".*\\b[D]?bNal") ||
                            shortcut.equals("Y")  ||  shortcut.matches(".*\\b[D]?(Me)?Tyr"))
                      commands.addElement("c " + rasmol_colors[2]);    // blue
                   else if (shortcut.equals("G")  ||  shortcut.matches(".*\\bGly"))
                      commands.addElement("c " + rasmol_colors[3]);    // light grey
                   else if (shortcut.equals("A")  ||  shortcut.matches(".*\\b[D]?(Me)?Ala"))
                      commands.addElement("c " + rasmol_colors[4]);    // dark grey
                   else if (shortcut.equals("H")  ||  shortcut.matches(".*\\b[D]?(Me)?His"))
                      commands.addElement("c " + rasmol_colors[5]);    // pale blue
                   else if (shortcut.equals("C")  ||  shortcut.matches(".*\\b[D]?(Me)?Cys")  ||
                            shortcut.equals("M")  ||  shortcut.matches(".*\\b[D]?(Me)?Met"))
                      commands.addElement("c " + rasmol_colors[6]);    // yellow
                   else if (shortcut.equals("S")  ||  shortcut.matches(".*\\b[D]?(Me)?Ser")  ||
                            shortcut.equals("T")  ||  shortcut.matches(".*\\b[D]?(Me)?Thr"))
                      commands.addElement("c " + rasmol_colors[7]);    // orange
                   else if (shortcut.equals("N")  ||  shortcut.matches(".*\\b[D]?(Me)?Asn")  ||
                            shortcut.equals("Q")  ||  shortcut.matches(".*\\b[D]?(Me)?Gln"))
                      commands.addElement("c " + rasmol_colors[8]);    // cyan
                   else if (shortcut.equals("L")  ||  shortcut.matches(".*\\b[D]?(Me)?Leu")  ||
                            shortcut.equals("I")  ||  shortcut.matches(".*\\b[D]?(Me)?Ile")  ||
                                                      shortcut.matches(".*\\b[D]?(Me)?Nle")  ||
                                                      shortcut.matches(".*\\b[D]?(Me)?Cha")  ||
                            shortcut.equals("V")  ||  shortcut.matches(".*\\b[D]?(Me)?Val")  ||
                                                      shortcut.matches(".*\\b[D]?(Me)?Nva"))
                      commands.addElement("c " + rasmol_colors[9]);    // green

                   else if (shortcut.equals("W")  ||  shortcut.matches(".*\\b[D]?(Me)?Trp"))
                      commands.addElement("c " + rasmol_colors[10]);    // pink
                   else if (shortcut.equals("P")  ||  shortcut.matches(".*\\b[D]?(Me)?Pro"))
                      commands.addElement("c " + rasmol_colors[11]);    // flesh
                   else
                      commands.addElement("c " + colortable[0]);
               }
               else if (a.symbol.equals("R")  &&  0 != (flags & USE_CLUSTAL))    /* shortcuts colored the CLUSTAL way */
               {
                   String shortcut = a.getStringProperty("A");
                   if (shortcut == null) 
                      commands.addElement("c " + colortable[0]);
                   else if (shortcut.equals("G")  ||  shortcut.matches(".*\\bGly"))
                      commands.addElement("c " + clustal_colors[0]);    // ORANGE
                   else if (shortcut.equals("P")  ||  shortcut.matches(".*\\b[D]?(Me)?Pro"))
                      commands.addElement("c " + clustal_colors[1]);    // YELLOW
                   else if (shortcut.equals("A")  ||  shortcut.matches(".*\\b[D]?(Me)?Ala")  ||
                            shortcut.equals("I")  ||  shortcut.matches(".*\\b[D]?(Me)?Ile")  ||
                            shortcut.equals("L")  ||  shortcut.matches(".*\\b[D]?(Me)?Leu")  ||
                            shortcut.equals("M")  ||  shortcut.matches(".*\\b[D]?(Me)?Met")  ||
                            shortcut.equals("F")  ||  shortcut.matches(".*\\b[D]?(Me)?Phe")  ||
                            shortcut.equals("W")  ||  shortcut.matches(".*\\b[D]?(Me)?Trp")  ||
                            shortcut.equals("V")  ||  shortcut.matches(".*\\b[D]?(Me)?Val"))
                      commands.addElement("c " + clustal_colors[2]);    // BLUE
                   else if (shortcut.equals("R")  ||  shortcut.matches(".*\\b[D]?(Me)?Arg")  ||
                            shortcut.equals("K")  ||  shortcut.matches(".*\\b[D]?(Me)?Lys"))
                      commands.addElement("c " + clustal_colors[3]);    // RED
                   else if (shortcut.equals("N")  ||  shortcut.matches(".*\\b[D]?(Me)?Asn")  ||
                            shortcut.equals("Q")  ||  shortcut.matches(".*\\b[D]?(Me)?Gln")  ||
                            shortcut.equals("S")  ||  shortcut.matches(".*\\b[D]?(Me)?Ser")  ||
                            shortcut.equals("T")  ||  shortcut.matches(".*\\b[D]?(Me)?Thr"))
                      commands.addElement("c " + clustal_colors[4]);    // GREEN
                   else if (shortcut.equals("C")  ||  shortcut.matches(".*\\b[D]?(Me)?Cys"))
                      commands.addElement("c " + clustal_colors[5]);    // PINK
                   else if (shortcut.equals("E")  ||  shortcut.matches(".*\\b[D]?(Me)?Glu")  ||
                            shortcut.equals("D")  ||  shortcut.matches(".*\\b[D]?(Me)?Asp"))
                      commands.addElement("c " + clustal_colors[6]);    // MAGENTA
                   else if (shortcut.equals("H")  ||  shortcut.matches(".*\\b[D]?(Me)?His")  ||
                            shortcut.equals("Y")  ||  shortcut.matches(".*\\b[D]?(Me)?Tyr"))
                      commands.addElement("c " + clustal_colors[7]);    // CYAN
                   else
                      commands.addElement("c " + colortable[0]);
               }
               else if (!a.symbol.equals("C"))     /* dark cyan */
               {
                  commands.addElement("c " + colortable[6]);
               }
               else
                  color_changed = false;
            }
            else if ((flags & USE_ATOM_COLORS) != 0  &&  (colors != null))
            {
               if (!colortable[colors[i]].equals(colortable[0]))
                  color_changed = true;
               else
                  color_changed = false;
               commands.addElement("c " + colortable[colors[i]]);
            }

            if (labels   != null  &&
               labels[i] != null  &&
               !labels[i].equals(""))        /* there is a label */
               symbolp = labels[i];
            else if (a.isotope != Atom.NATURAL  &&  a.symbol.equals("H"))
            {
                if (a.isotope == 2)      symbolp = "D";
                else if (a.isotope == 3) symbolp = "T";
                else                     symbolp = atom_symbol;
                System.err.println("Setting atom symbol = " + symbolp);
            }
            else
               symbolp = atom_symbol;

            if (charge    != 0                ||        /* tricky case */
                a.radical != Atom.NO_RADICAL  ||
                (a.isotope != Atom.NATURAL  &&  !a.symbol.equals("H")))
            {
               commands.addElement("t " + direction + "01 " +
                                 worldToMetaX(a.x) + " " +
                                 worldToMetaY(a.y) + " " +
                                 symbolp + " " +
                                 AttributesToString(charge,
                                                    a.radical,
                                                    a.isotope));
            }
            else  /* ordinary case */
               commands.addElement("t " + direction + "0 " +
                                 worldToMetaX(a.x) + " " +
                                 worldToMetaY(a.y) + " " +
                                 symbolp);
            if (color_changed  &&  (flags & USE_ATOM_COLORS) != 0)
            {
               commands.addElement("c " + colortable[0]);
            }
         }
      }

      /* draw stereo uncertainty of double bonds */
      for (int i=0; i<m.bonds.length; i++)
      {
         b = m.bonds[i];
         if (b.stereo_symbol == Bond.DBL_EITHER ||
             (b.bond_type == Bond.SINGLE  &&  b.stereo_symbol == Bond.EITHER))
         {
            x1 = worldToMetaX(m.atoms[b.atoms[0]-1].x);
            y1 = worldToMetaY(m.atoms[b.atoms[0]-1].y);
            x2 = worldToMetaX(m.atoms[b.atoms[1]-1].x);
            y2 = worldToMetaY(m.atoms[b.atoms[1]-1].y);
            commands.addElement("t 0 " + ((x1+x2)/2) + " " +
                                         ((y1+y2)/2) + " ?");
         }
      }

      String[] result = new String[commands.size()];
      for (int i=0; i<commands.size(); i++)
         result[i] = commands.elementAt(i).toString();

      return (result);
   }
}  // DepictMolecule
