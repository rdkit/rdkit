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

import java.awt.BasicStroke;
import java.awt.Stroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;
import java.util.Vector;

public class MetaDraw {
	static Map<String, Font> fcache = new HashMap<String, Font>();
	static Object cache_lock = new Object();

	// define defaults
	static double defaultMinFontSize = 5.0;
	static float defaultStrokeWidth = 1.0f;

	/**
	 * Returns preferred plot size of the drawing given by plotCommands[].
	 */
   public static Dimension preferredPlotSize(String[] plotCommands)
   {
		Dimension result = new Dimension(100, 100);

      if (plotCommands == null) return (result);

		int xExt, yExt;

		int dpi = 72;
      try { Toolkit.getDefaultToolkit().getScreenResolution(); }
      catch (Exception e) { /* catch HeadLessException if JDK 1.4 */ }

      for (int i=0; i<plotCommands.length; i++)
      {
			StringTokenizer st = new StringTokenizer(plotCommands[i]);
			String opcode = st.nextToken();
			if (opcode.equals("e")) // size of world coordinate window
			{
				xExt = (int) (Integer.parseInt(st.nextToken()) * (dpi / (72.0 * 20)));
				yExt = (int) (Integer.parseInt(st.nextToken()) * (dpi / (72.0 * 20)));
				result = new Dimension(Math.max(100, xExt), Math.max(100, yExt));
			}
		}

		return (result);
	}

	/**
	 * Parse <code>symbol</code> into element segments.
	 * 
     * For example, "C2H5OH" is parsed as "C2" | "H5" | "O" | "H".
     * This method is used to implement the drawing path of chemical
     * symbol labels.
	 */
    private static String[] parseSymbol(String symbol)
    {
		Vector<String> segments = new Vector<String>();
		int index = 0;
		StringBuffer tmp = new StringBuffer();
        while (index < symbol.length())
        {
			// grab first character
			tmp.append(symbol.charAt(index++));
			// grab any lower case subsequent characters
            while (index < symbol.length()  &&
                   Character.isLowerCase(symbol.charAt(index)))
            {
				tmp.append(symbol.charAt(index++));
			}
			// grab any subsequent digits
            while (index < symbol.length()  &&
                   Character.isDigit(symbol.charAt(index)))
            {
				tmp.append(symbol.charAt(index++));
			}

			segments.add(tmp.toString());
			tmp.setLength(0);
		}

		String result[] = new String[segments.size()];
		for (int i = 0; i < result.length; i++)
			result[i] = (String) segments.get(i);
		return result;
	}

   /**
    * Draw the text defined by the tokens in <code>st</code> to the
    * Graphics context g.
    * 
    * The other parameters define the coordinate transformation to be used.
    */
   private static void drawText(Graphics g, StringTokenizer st,
                                int xOffset, int yOffset,
                                int xOrg, int yOrg,
                                double scale)
   {
      Graphics2D g2 = null;
      if (g instanceof Graphics2D) g2 = (Graphics2D)g;

      // text symbols consist of a center piece
      // plus a number of surrounding strings, e.g. charge or mass
      String mode = st.nextToken().toLowerCase();

      FontMetrics fm = g.getFontMetrics();
      // set default width of first character if there is no center symbol
      float fw = fm.charWidth('N');
      float fh = fm.getAscent();

      // Coordinate where the center of the first character of the symbol
      // shall be placed.
      double x1 = Integer.parseInt(st.nextToken());
      double y1 = Integer.parseInt(st.nextToken());

      float xoff = xOffset;
      float yoff = yOffset;
      double xDraw, yDraw;

      char direction = 'e'; // default: write eastward, i.e. left-to-right
      if (Character.isLowerCase(mode.charAt(0)))
      {
         direction = mode.charAt(0);
         mode = mode.substring(1);
      }
      double symbolHeight = fh;
      double symbolWidth = fw;

      // Depending on the mode, there may be more than one string.
      // each digit defines the position of a string that is expected
      // in the token stream.
      if (mode.substring(0, 1).equals("0")) // default string
      {
         // String to be used for atomic symbol
         String symbol = st.nextToken();
         if (symbol == null  ||  symbol.trim().equals("")) return;
         // System.err.println("drawing '" + symbol + "' at center");
         // Parse it into pieces
         String symbolParts[] = parseSymbol(symbol);
         fw = fm.stringWidth(symbol.substring(0, 1));
         Font font = fm.getFont();
         if (g2 != null)
         {
            java.awt.font.GlyphVector gv = font.createGlyphVector(g2.getFontRenderContext(), symbol.substring(0,1));
            Shape shape = gv.getGlyphVisualBounds(0);
            fh = (float) shape.getBounds2D().getHeight();
            fw = (float) shape.getBounds2D().getWidth();
         }
         symbolWidth = fm.stringWidth(symbol);
         // Center on first character
         xoff -= fw/2.0; yoff += fh/2.0;

         // The default string is drawn with the reference point at the
         // center of the first character of the first symbol part. This
         // is not 100% OK for drawing more than one character parts
         // from right to left but it works.
         //
         // The symbol parts are printed left to right each.
         // The print direction is implemented by calculating offsets
         // before and after printing the string.
         double shiftWidth = 0;
         double shiftHeight = 0;
         for (int k=0; k<symbolParts.length; k++)
         {
            float sw = 0, sh = 0;
            symbol = symbolParts[k];
            float xDirOffset = 0;
            float yDirOffset = 0;
            if (k == 0)
            {
               symbolWidth = fm.stringWidth(symbol);
               if (direction == 'w') shiftWidth = (float) fm.stringWidth(symbol);
            }
            switch (direction)
            {
               case 'e':
                  xDirOffset = 0;
                  yDirOffset = 0;
                  break;
               case 's':
                  xDirOffset = 0;
                  break;
               case 'w':
                  xDirOffset = -(float) fm.stringWidth(symbol);
                  yDirOffset = 0;
                  break;
               case 'n':
                  xDirOffset = 0;
                  yDirOffset = 0;
                  break;
            }
            // reference for text in corners is only first symbol part
            // if string is not printed left-to-right.
            for (int j=0; j<symbol.length(); j++)
            {
               String s1 = symbol.substring(j, j + 1);
               xDraw = shiftWidth+xDirOffset +
               xoff+((x1-xOrg)*scale)+sw;
               yDraw = shiftHeight+yDirOffset+
               yoff+((y1-yOrg)*scale)+sh;
               // numbers are automatically subscripted
               if (Character.isDigit(s1.charAt(0))) yDraw += fh/2;
               if (g2 != null)
                  g2.drawString(s1, (float) xDraw, (float) yDraw);
               else
                  g.drawString(s1, (int)(xDraw+0.5), (int)(yDraw+0.5));
               sw += fm.stringWidth(s1);
            }
            switch (direction)
            {
               case 'e':
                  shiftWidth += fm.stringWidth(symbol) + 0.1 * fw;
                  break;
               case 's':
                  shiftHeight += 1.1 * fh;
                  break;
               case 'w':
                  shiftWidth -= fm.stringWidth(symbol) + 0.1 * fw;
                  break;
               case 'n':
                  shiftHeight -= 1.1 * fh;
                  break;
            }
         }
      }

      // draw surrounding strings
      while (!(mode = mode.substring(1)).equals(""))
      {
         String ss = st.nextToken();

         xDraw = xoff + ((x1 - xOrg) * scale);
         yDraw = yoff + ((y1 - yOrg) * scale);
         if (mode.charAt(0) == '1') // upper right
         {
            xDraw += symbolWidth;
            yDraw -= fh;
         }
         else if (mode.charAt(0) == '2')  // lower right
         {
            xDraw += symbolWidth;
            yDraw += fh / 2;
         }
         else if (mode.charAt(0) == '3')  // lower left
         {
            xDraw -= fw;
            yDraw -= fh;
         }
         else if (mode.charAt(0) == '4')  // upper left
         {
            xDraw -= fw;
            yDraw += fh / 2;
         }

         if (g2 != null)
         {
            g2.drawString(ss, (float) xDraw, (float) yDraw);
         }
         else
         {
            g.drawString(ss, (int) xDraw, (int) yDraw);
         }
      }

   }

   public static void depictCommands(Graphics g,
                                     String[] plotCommands,
                                     Rectangle aRect)
   {
      depictCommands(g, plotCommands, aRect, true, false);
   }

   // Delegates processing to more general method with suitable defaults
   public static void depictCommands(Graphics g,
                                     String[] plotCommands,
                                     Rectangle aRect,
                                     boolean colored,
                                     boolean forColorBlind)
   {
      depictCommands(g, plotCommands, aRect, colored, forColorBlind, defaultMinFontSize, defaultStrokeWidth);
   }

   public static void depictCommands(Graphics g,
                                     String[] plotCommands,
                                     Rectangle aRect,
                                     boolean colored,
                                     boolean forColorBlind,
                                     double minFontSize,
                                     float strokeWidth)
   {
      int xExt, yExt, xOrg, yOrg;
      int xoff, yoff;
      double width, height;
      double scale;

      int dpi = 72;
      try { Toolkit.getDefaultToolkit().getScreenResolution(); }
      catch (Exception e) { /* catch HeadLessException if JDK 1.4 */ }

      Graphics2D g2 = null;
      BasicStroke defaultStroke = new BasicStroke(strokeWidth);
      BasicStroke boldStroke = new BasicStroke(strokeWidth*2.0f);
      if (g instanceof Graphics2D)
      {
         g2 = (Graphics2D) g;
         g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                             RenderingHints.VALUE_ANTIALIAS_ON);
         g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS,
                             RenderingHints.VALUE_FRACTIONALMETRICS_ON);
         g2.setStroke(defaultStroke);
      }
      width = aRect.width;
      height = aRect.height;
      xoff = aRect.x;
      yoff = aRect.y;
      xExt = yExt = 100;
      xOrg = yOrg = -50;

      scale = width / xExt;
      if (scale > height/yExt) scale = height/yExt;


      Line2D.Double line = new Line2D.Double(); // out of loop for performance
      GeneralPath path = new GeneralPath();
      int currentColor = 0;
      if (plotCommands == null  ||  plotCommands.length == 0) return;
      for (int i=0; i<plotCommands.length; i++)
      {
         String opcode;

         try 
         {
            StringTokenizer st = new StringTokenizer(plotCommands[i]);
            opcode = st.nextToken();
            if (opcode.equals("e")) // size of world coordinate window
            {
               xExt = Integer.parseInt(st.nextToken());
               yExt = Integer.parseInt(st.nextToken());
               scale = width / xExt;
               if (scale > height/yExt) scale = height/yExt;
               if (scale > dpi/(72.0*20)) scale = dpi/(72.0*20); // not too big
               xoff = xoff + (int) (width - scale * xExt) / 2;
               yoff = yoff + (int) (height - scale * yExt) / 2;
            }
            else if (opcode.equals("o"))  // offset of world origin
            {
               xOrg = Integer.parseInt(st.nextToken());
               yOrg = Integer.parseInt(st.nextToken());
            }
            else if (opcode.equals("c"))  // set a color of if mode is colored
            {
               if (colored)
               {
                  int col = Integer.parseInt(st.nextToken().substring(2),16);
                  currentColor = col;
                  g.setColor(new Color(col));
               }
            }
            else if (opcode.equals("f"))  // set a font
            {
               String fname = st.nextToken();
               int size = Integer.parseInt(st.nextToken()); // in TWIPS
               int fontStyle = Font.BOLD;
               // try to parse entry for font style
               try
               {
                  String style = st.nextToken();
                  if (style != null)
                  {
                     if (style.equalsIgnoreCase("p")) { fontStyle = Font.PLAIN; }
                     else if (style.equalsIgnoreCase("i"))
                     {
                        fontStyle = Font.ITALIC;
                     }
                     else if (style.equalsIgnoreCase("bi"))
                     {
                        fontStyle |= Font.ITALIC;
                     }
                  }
               }
               catch (NoSuchElementException e)
               {
                  // just do nothing
               }
               double dsize = (double) size / 20; // in points
               dsize = (dsize*(scale/(dpi/(72.0*20))));     // to fit box scale
               // if (dsize < 7) dsize = 7; // ensure visibility
               if (dsize < minFontSize) dsize = minFontSize;

               Font f;
               String key = (int) (dsize * 3 + 0.5) + "x" + (fontStyle + 1);
               synchronized (cache_lock)
               {
                  if (fcache.containsKey(key))
                  {
                     f = fcache.get(key);
                     g.setFont(f);
                  }
                  else
                  {
                     Font ftmp = new Font(fname, fontStyle, 5);
                     f = ftmp.deriveFont((float) dsize);
                     fcache.put(key, f);
                     g.setFont(f);
                  }
               }
            }
            else if (opcode.equals("l"))  // draw a line
            {
               int x1 = Integer.parseInt(st.nextToken());
               int y1 = Integer.parseInt(st.nextToken());
               int x2 = Integer.parseInt(st.nextToken());
               int y2 = Integer.parseInt(st.nextToken());

               if (g2 != null)
               {
                  line.setLine(xoff+((x1-xOrg)*scale),
                  yoff+((y1-yOrg)*scale),
                  xoff+((x2-xOrg)*scale),
                  yoff+((y2-yOrg)*scale));
                  try
                  {
                     if (forColorBlind  &&  currentColor != 0)
                     {
                         Stroke stroke = g2.getStroke();
                         g2.setStroke(boldStroke);
                         g2.draw(line);
                         g2.setStroke(stroke);
                     }
                     else
                         g2.draw(line);
                  }
                  catch (Throwable t)
                  {
                     t.printStackTrace();
                     System.err.println("line = (" + line.getX1() + "," + line.getY1() + ")-(" +
                                                     line.getX1() + "," + line.getY2() + ")");
                  }
               }
               else
                  g.drawLine(xoff+(int)((x1-xOrg)*scale),
                             yoff+(int)((y1-yOrg)*scale),
                             xoff+(int)((x2-xOrg)*scale),
                             yoff+(int)((y2-yOrg)*scale));
            }
            // pf 3 -435 276 -269 633 -209 599
            else if (opcode.equals("pf")) // draw a filled polygon
            {
               int x, y;

               int n = Integer.parseInt(st.nextToken());
               int[] xa = new int[n];
               int[] ya = new int[n];
               double[] dxa = new double[n];
               double[] dya = new double[n];

               for (int j=0; j<n; j++)
               {
                  x = Integer.parseInt(st.nextToken());
                  y = Integer.parseInt(st.nextToken());
                  xa[j] = xoff + (int) ((x - xOrg) * scale);
                  ya[j] = yoff + (int) ((y - yOrg) * scale);
                  dxa[j] = xoff + ((x - xOrg) * scale);
                  dya[j] = yoff + ((y - yOrg) * scale);
               }

               if (g2 != null)
               {
                  path.reset();
                  for (int j=1; j<n; j++)
                  {
                     line.setLine(dxa[j - 1], dya[j - 1], dxa[j], dya[j]);
                     path.append(line, true);
                  }
                  line.setLine(dxa[n - 1], dya[n - 1], dxa[0], dya[0]);
                  path.closePath();
                  // g2.draw(path);
                  g2.fill(path);
               }
               else
               {
                  g.drawPolygon(xa, ya, n);
                  g.fillPolygon(xa, ya, n);
               }
            }
            else if (opcode.equals("t")) // draw a text symbol
            {
               // draw the text string defined by this line
               if (forColorBlind  &&  currentColor != 0)
               {
                  Font oldFont = g.getFont();
                  g.setFont(oldFont.deriveFont(Font.ITALIC));
                  drawText(g, st, xoff, yoff, xOrg, yOrg, scale);
                  g.setFont(oldFont);
               }
               else
                  drawText(g, st, xoff, yoff, xOrg, yOrg, scale);
            } 
            else if (opcode.equals("tl")) // draw a text symbol
            {
               // draw the text label defined by this line
               drawLabel(g, st, xoff, yoff, xOrg, yOrg, scale);
            } 
            else if (opcode.equals("ac")) // show atom coordinates
            {
               // Just for information => NOP
            }
            else
            {
               // unknown command => NOP
            }
         }
         catch (Exception e)
         {
            e.printStackTrace(System.err);
            System.err.println("problem '" + e + "' decoding command '" + plotCommands[i] +"'");
         }
      }
   }

	/*
	 * Method to draw a centered text label (added for SRNA depiction)
	 */

	private static void drawLabel(Graphics g, StringTokenizer st, int xOffset,
			int yOffset, int xOrg, int yOrg, double scale) {
		Graphics2D g2 = null;
		if (g instanceof Graphics2D)
			g2 = (Graphics2D) g;

		// text symbols consist of a center piece
		// plus a number of surrounding strings, e.g. charge or mass
		String mode = st.nextToken().toLowerCase();

		FontMetrics fm = g.getFontMetrics();
		// set default width of first character if there is no center symbol
		float fw = fm.charWidth('N');
		float fh = fm.getAscent();

		// Coordinate where the center of the first character of the symbol
		// shall be placed.
		double x1 = Integer.parseInt(st.nextToken());
		double y1 = Integer.parseInt(st.nextToken());

		float xoff = xOffset;
		float yoff = yOffset;
		double xDraw, yDraw;

		char direction = 'e'; // default: write eastward, i.e. left-to-right
		if (Character.isLowerCase(mode.charAt(0))) {
			direction = mode.charAt(0);
			// mode = mode.substring(1);
		}
		// double symbolHeight = fh;
		// double symbolWidth = fw;

		// String to be used for label
		String symbol = st.nextToken();
		if (symbol == null || symbol.trim().equals(""))
			return;
		fw = fm.stringWidth(symbol);
                java.awt.geom.Rectangle2D rect = fm.getStringBounds(symbol, g);
                fw = (float)rect.getWidth();
		Font font = fm.getFont();
		// Center on first character
		xoff -= fw / 2.0;
		yoff += fh / 2.0;
		xDraw = xoff + ((x1 - xOrg) * scale);
		yDraw = yoff + ((y1 - yOrg) * scale);
		if (g2 != null)
			g2.drawString(symbol, (float) xDraw, (float) yDraw);
		else
			g.drawString(symbol, (int) (xDraw + 0.5), (int) (yDraw + 0.5));

	}

}
