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
package avalon.jsp.servlets;

import org.apache.log4j.Logger;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.image.FilteredImageSource;
import java.awt.image.ImageFilter;
import java.awt.image.ImageProducer;
import java.awt.image.RGBImageFilter;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.font.FontRenderContext;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.security.MessageDigest;

import java.util.ArrayList;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Properties;

import javax.servlet.ServletException;
import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import avalon.tools.Base64;

import jni.JNIMatch;
import jni.JNIDepict;
import novartis.chemistry.molecule.MOLCanvas;
import novartis.chemistry.molecule.Molecule;
import novartis.chemistry.molecule.MoleculeDepicter;
import novartis.chemistry.molecule.StandardMOL;
import novartis.chemistry.molecule.SmilesParser;
import novartis.chemistry.molecule.Layout;
import novartis.utilities.MetaDraw;
import novartis.utilities.FortranInputStream;
import avalon.jni.JNISmi2Mol;

import java.util.Iterator;
import javax.imageio.ImageWriter;
import javax.imageio.ImageIO;

/**
 * This servlet is used to convert MOL file input or SMILES parameters to
 * an image format.
 * <br>
 * The idea is to separate textual content generation from image rendering,
 * since image rendering may require a richer JVM to process it. 
 */
public class Mol2ImageServlet extends HttpServlet
{
	static String buildNumber = "DEFAULT";
	static {
		try {
		//load a properties file from class path, inside static method
		Properties prop = new Properties();
		InputStream ps = Mol2ImageServlet.class.getClassLoader().getResourceAsStream("build.properties");
		if (ps != null) {
			prop.load(ps);
			buildNumber = prop.getProperty("build.number");
			ps.close();
		}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	/* Get actual class name to be printed on */
    static Logger log = Logger.getLogger(Mol2ImageServlet.class.getName());
    
    public void doGet(HttpServletRequest req, HttpServletResponse res)
        throws ServletException, IOException
    {
        doService(req, res, false);   
    }
    public void doPost(HttpServletRequest req, HttpServletResponse res)
        throws ServletException, IOException
    {
        doService(req, res, true);   
    }

    // create an LRU cache with initial capacity 100 and load factor 0.8 ordered by last accessed
    static MOLCache<String, byte[]> cache = new MOLCache<String, byte[]>(100, 0.8f, true);

    static int nPNG = 0;
    static int nMD5 = 0;
    static int nGIF = 0;
    static int nWithScaffold = 0;
    static int nWithTplMOL = 0;
    static int nFPLabels = 0;
    static int nHEX = 0;
    static int nBIN = 0;
    static int nLST = 0;
    static int nPLT = 0;
    static int nFlags = 0;
    static int nB64 = 0;
    static int nMOL = 0;

    public void init()
    {
        ImageIO.setUseCache(false);
    }

    static int nHealth = 0;

    public synchronized void reportHealth(HttpServletResponse response)
        throws IOException
    {
        nHealth++;
        response.setContentType("text/xml");
        PrintWriter out = response.getWriter();
        out.println("<?xml version='1.0' encoding='UTF-8' standalone='yes'?>");
        Runtime rt = Runtime.getRuntime();
        out.println("<health>");
        out.println("    <heapSize>" + rt.totalMemory() + "</heapSize>");
        out.println("    <heapMaxSize>" + rt.maxMemory() + "</heapMaxSize>");
        out.println("    <heapFreeSize>" + rt.freeMemory() + "</heapFreeSize>");
        rt.gc();
        out.println("    <afterGC total='" + rt.totalMemory() + "' free='" + rt.freeMemory() + "' used='" + (rt.totalMemory()-rt.freeMemory()) + "'/>");
        out.println("    <molCache size='" + cache.size() + "'/>");
        out.println("    <counters");
        out.println("       nHealth='" + nHealth + "' ");
        out.println("       nPNG='" + nPNG + "' ");
        out.println("       nMD5='" + nMD5 + "' ");
        out.println("       nGIF='" + nGIF + "' ");
        out.println("       nBIN='"  + nBIN + "' ");
        out.println("       nHEX='"  + nHEX + "' ");
        out.println("       nLST='"  + nLST + "' ");
        out.println("       nPLT='" + nPLT + "' ");
        out.println("       nFlags='" + nFlags + "' ");
        out.println("       nMOL='" + nMOL + "' ");
        out.println("       nB64='" + nB64 + "' ");
        out.println("       nFPLabels='" + nFPLabels + "' ");
        out.println("       nWithTplMOL='" + nWithTplMOL + "' ");
        out.println("       nWithScaffold='" + nWithScaffold + "' " + "/>");
        out.println("    <version>1.2.0</version>");
        out.println("    <build>"+buildNumber+"</build>");
        out.println("</health>");
        out.flush();
    }

    /**
     * Handles a GET and POST requests by rendering SMILES (with optional coordinates) or MOL files
     * and sending the picture to the caller.
     *
     * Note: The rendering does not reqire a session ID because it does not
     * connect to data but just transforms it.
     *
     * @param request an HttpServletRequest
     * @param response an HttpServletResponse
     * @exception IOException
     */
    public void doService(HttpServletRequest request, HttpServletResponse response, boolean isPostRequest) 
        throws IOException
    {
        // Get the pathInfo. The extension determines the result format:
        // E.g., .gif yields a GIF file and .png a Publig Network Graphics image.
        String path = request.getPathInfo(); 
        if (path.endsWith("health.xml"))
        {
            reportHealth(response);
            return;
        }

        String footer = request.getParameter("footer");       // Footer string to be rendered in molecule font bottom/centered
        String title = request.getParameter("title");         // Title string to be rendered in larger font top/centered

        // collect all possible parameter values from request (even if not needed for certain code paths)
        String smiles = request.getParameter("smiles");         // SMILES to be rendered
        String reactionSmiles = null;
        if (smiles != null)
        {
            smiles = smiles.trim();
            // remove valence definition if any
            smiles = smiles.replaceAll("\\;v\\d+", "");
            // convert Rn syntax of R-groups to mass notation
            smiles = smiles.replaceAll("\\[R(\\d+)\\]", "[$1R]");
            // fix odd R-group syntax
            smiles = smiles.replaceAll("\\[([A-Z][a-z]*)\\*\\]", "[$1](*)");
            log.debug("pruned smiles = '" + smiles + "'");
        }
        if (smiles != null  &&  smiles.trim().matches("[^>]+>[^>]*>[^>]+")) reactionSmiles = smiles;
        boolean useChains = "true".equalsIgnoreCase(request.getParameter("chains"));  // abbreviate SMILES chains
        String molfile = request.getParameter("molfile");       // get MOL file if any
        String scaffold = request.getParameter("scaffold");     // SMILES to be used for consistent orientation
        if (scaffold != null  &&  scaffold.trim().length() == 0) scaffold = null;  // completely clear empty templates
        if (scaffold != null) {nWithScaffold++; log.debug("scaffold =\n" + scaffold);}
        String tplMOL = request.getParameter("tplMOL");         // the scaffold could also be defined as a MOL file parameter
        if (tplMOL != null  &&  tplMOL.trim().length() == 0) tplMOL = null;  // completely clear empty templates
        if (tplMOL != null) {
            nWithTplMOL++;
        }
        // dimensions of image if any
        String swidth = request.getParameter("width");          // width of image box
        if (swidth == null) swidth = request.getParameter("w"); // use shorter version if not yet defined
        int width = 300;
        if (swidth != null) try {width = Integer.parseInt(swidth);} catch (RuntimeException e) {/* NOP */}
        String sheight = request.getParameter("height");         // height of image box
        if (sheight == null) sheight = request.getParameter("h");
        int height = 250;
        if (sheight != null) try {height = Integer.parseInt(sheight);} catch (RuntimeException e) {/* NOP */}

        String md5 = request.getParameter("md5");               // return a picture as indexed by the MD5 of the MOL file posted previously

        String nbitsString = request.getParameter("nbits");     // return a hex-coded bit pattern of nbits bits.
        int nbits = 512;                                        // Note nbits == 2^N. 512 is a reasonable value yielding an average density of ~15%
        try {nbits = Integer.parseInt(nbitsString);} catch (Exception e) {/* NOP */} 

        boolean asQuery = false;
        if ("true".equalsIgnoreCase(request.getParameter("asquery"))) asQuery = true;

        // collect fingerprint coloring info
        boolean useFPLabels = false;
        if ("true".equalsIgnoreCase(request.getParameter("usefplabels"))) { nFPLabels++; useFPLabels = true;}
        String bits = request.getParameter("bits");
        byte[] tmpfp = new byte[nbits/8];
        byte[] firstBitFP = new byte[nbits/8];
        int[] bitList = null;
        boolean forColorBlind = "true".equalsIgnoreCase(request.getParameter("cb"));
        if (bits != null)
        {
            String[] bitStringList = bits.split(",");
            bitList = new int[bitStringList.length];
            int ii = 0;
            for (int i=0; i<bitStringList.length; i++)
            {
                try
                {
                    int bit = Integer.parseInt(bitStringList[i]);
                    if (bit >= nbits) continue;  // skip bits that are beyond FP range
                    bitList[ii] = bit;
                    tmpfp[bit/8] |= (byte)(1<<(bit%8));
                    if (ii == 0)
                    {
                        firstBitFP[bit/8] |= (byte)(1<<(bit%8));
                    }
                    ii++;
                }
                catch (Exception e)
                {
                    /* NOP */
                }
            }
            if (ii < bitStringList.length)  // strip entries that aren't fingerprints
            {
                int[] newBitList = new int[ii];
                for (int iii=0; iii<ii; iii++) newBitList[iii] = bitList[iii];
                bitList = newBitList;
            }
        }

        // collect a few more parameters
        String colored = request.getParameter("colored");
        String tpl_coords = request.getParameter("tpl_coords");
        boolean hasBorder = null != request.getParameter("b");
        String tmp = request.getParameter("bg");
        String transparent = request.getParameter("t");
        String flagString = request.getParameter("flags");
        if (flagString != null  &&  flagString.matches(".*\\bchains\\b.*")) useChains = true;
        int flags = 0;
        if (flagString != null  &&  flagString.matches(".*\\baminoacids\\b.*")) flags |= JNISmi2Mol.APPLY_AMINO_ACIDS;
        if (flagString != null  &&  flagString.matches(".*\\bprotecting-groups\\b.*")) flags |= JNISmi2Mol.APPLY_PROTECTING_GROUPS;
        if (flagString != null  &&  flagString.matches(".*\\bextended\\b.*")) flags |= JNISmi2Mol.EXTENDED_SHORTCUTS;
        if (flagString != null  &&  flagString.matches(".*\\bnon-standard\\b.*")) flags |= JNISmi2Mol.NON_STANDARD_SHORTCUTS;
        if (flagString != null  &&  flagString.matches(".*\\bcatch-all\\b.*")) flags |= JNISmi2Mol.CATCH_ALL_SHORTCUTS;
        if (flagString != null  &&  flagString.trim().length() > 0) log.debug("flagString = '" + flagString + "', flags = " + flags);
        if (flags != 0) nFlags++;

        long now = System.currentTimeMillis();
        response.setDateHeader("Expires", now + 24*60*60*1000);

        // normal processing
        // Prepare template for coloring and orientation
        if (tplMOL == null  &&  scaffold != null  &&  !scaffold.trim().equals(""))
        {
            synchronized (Mol2ImageServlet.class)
            {
                tplMOL = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(scaffold, JNISmi2Mol.EXPECT_SMARTS);
                if (tplMOL == null  || tplMOL.trim().length() == 0) log.info("Could not convert scaffold '" + scaffold + "' to tplMOL");
            }
        }
        // Convert SMILES to MOL file if any
        if (molfile == null  &&  smiles != null  && !smiles.equals(""))
        {
            if (useChains)
            {
                smiles = abbreviateChains(smiles);
                log.debug("after useChains: smiles = '" + smiles + "'");
            }
            synchronized (Mol2ImageServlet.class)
            {
                if (tplMOL != null  &&  !tplMOL.trim().equals(""))
                {
                    if (tpl_coords != null  &&  tpl_coords.equalsIgnoreCase("rotate"))
                        molfile = JNIDepict.getDepictor().smiToMOLWithTemplate(smiles, tplMOL);
                    else if (tpl_coords != null  &&  tpl_coords.equalsIgnoreCase("force"))
                        molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithTemplate(smiles, tplMOL);
                    else
                        molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(smiles, flags);
                }
                else
                {
                    molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(smiles, flags);
                }
            }
        }
        if (molfile == null) {
            molfile = "";
            log.info("smiles = " + smiles);
            log.info("tpl_coords = " + tpl_coords);
            log.info("flags = " + flags);
            log.info("tplMOL = " + tplMOL);
        }
        if (path.endsWith(".smi"))
        {
            response.setContentType("text/plain");
            PrintWriter out = response.getWriter();
            synchronized (Mol2ImageServlet.class)
            {
                if (molfile.length() < 30)
                    out.println("Broken MOL file\n" + molfile);
                else
                    out.println(JNISmi2Mol.getSmi2Mol().getSMILESFromCT(molfile));
            }
            return;
        }

        MOLCanvas mol = null;
        String[] labels = null;
        if (molfile.length() > 0)
        {
            mol = new MOLCanvas(molfile.toString(), width, height, false);
            if (tplMOL != null  &&  ("yes".equalsIgnoreCase(colored)  ||  "true".equalsIgnoreCase(colored)))     // recompute depiction commands to color match
            {
                try
                {
                    BufferedReader br = new BufferedReader(new StringReader(molfile.toString()));
                    FortranInputStream in = new FortranInputStream(br);
                    StandardMOL m = new StandardMOL();
                    m.readMOLFile(in);
                    // color matching atoms in blue
                    JNIMatch matcher = JNIMatch.getMatcher();
                    int[] colorIndices = matcher.getMatchColoring(molfile.toString(), tplMOL, 2);
                    mol.setMOLFile(molfile.toString(), colorIndices);
                }
                catch (Throwable e)
                {
                    System.err.println("Exception " + e + " ignored");
                }
            }
            else if (nbits > 0  &&  bitList != null  &&  smiles != null)
            {
                int[] colorIndices = new int[smiles.length()];
                JNIMatch matcher = JNIMatch.getMatcher();
                int mask = 0;
                synchronized (Mol2ImageServlet.class)
                {
                    for (int t=0; t<31; t++)
                    {
                        // byte[] molfp = new byte[tmpfp.length];
                        byte[] molfp = 
                        matcher.getFingerprintFromCT(molfile, false, (1<<t), tmpfp.length);
                        for (int i=0; i<tmpfp.length; i++)
                            if (tmpfp[i] != 0  &&  0 != (molfp[i]&tmpfp[i]))
                            {
                                mask |= (1<<t);
                            }
                    }
                }

                HashSet<String> fpClassSet = new HashSet<String>();
                ArrayList<String> fpBitList = new ArrayList<String>();
                // grab colors of first (least likely) bit first
                for (int t=0; t<31; t++)
                {
                    int[] tmpColors = new int[colorIndices.length];
                    if (0 == (mask&(1<<t))) continue;
                    synchronized (Mol2ImageServlet.class)
                    {
                        matcher.setFingerprintColorsForCT(molfile, firstBitFP, false, 1<<t, tmpColors);
                    }
                    // use red color for all changes in priority bit
                    for (int i=0; i<tmpColors.length; i++)
                        if (tmpColors[i] > 0  &&  colorIndices[i] == 0)
                        {
                            for (int ibit=0; ibit<firstBitFP.length*8; ibit++)
                                if ((firstBitFP[ibit/8]&(byte)(1<<(ibit%8))) != 0)
                                {
                                    String bitString = JNIMatch.maskToString(1<<t) + "(" + ibit + ")";
                                    if (!fpClassSet.contains(bitString))
                                    {
                                        fpClassSet.add(bitString);
                                        fpBitList.add(bitString);
                                    }
                                }
                            colorIndices[i] = 1;        // red
                        }
                }

                for (int i=0; i<tmpfp.length; i++) tmpfp[i] ^= firstBitFP[i];   // clear away first bit
                // now, get additional colors for other bits
                for (int t=0; t<31; t++)
                {
                    int[] tmpColors = new int[colorIndices.length];
                    if (0 == (mask&(1<<t))) continue;
                    synchronized (Mol2ImageServlet.class)
                    {
                        matcher.setFingerprintColorsForCT(molfile, tmpfp, false, 1<<t, tmpColors);
                    }
                    // use red color for all forced atoms and blue for changed atoms
                    for (int i=0; i<tmpColors.length; i++)
                    {
                        if (tmpColors[i] > 0)
                            for (int ibit=0; ibit<tmpfp.length*8; ibit++)
                                if ((tmpfp[ibit/8]&(byte)(1<<(ibit%8))) != 0)
                                {
                                    String bitString = JNIMatch.maskToShortString(1<<t);
                                    if (!fpClassSet.contains(bitString))
                                    {
                                        fpClassSet.add(bitString);
                                        fpBitList.add(bitString);
                                    }
                                }
                        if (tmpColors[i] == 1  &&  colorIndices[i] != 1)
                        {
                            colorIndices[i] = 2;        // blue
                        }
                        else if (tmpColors[i] == 2  &&  colorIndices[i] == 0)   // changed atom
                        {
                            colorIndices[i] = 7;        // dark yellow
                        }
                    }
                }
                mol.setMOLFile(molfile.toString(), colorIndices);
                if (useFPLabels) labels = (String[])fpBitList.toArray(new String[0]);
            }
        }
        else
        {
            mol = new MOLCanvas(null, width, height, false);
        }
        mol.setBackground(java.awt.Color.white);
        mol.setOpaque(true);
        mol.setSize(width,height);

        Color backColor = Color.white;
        if (tmp != null)
            try {backColor = new Color(Integer.parseInt(tmp,16));}
            catch (Exception e) //{/* Ignore */}
            {
                System.err.println("bg = " + tmp);
                e.printStackTrace();
            }

        if (path.endsWith(".gif")) nGIF++;
        else if (path.endsWith(".png")) nPNG++;
        else if (path.endsWith(".hex")) nHEX++;
        else if (path.endsWith(".bin")) nBIN++;
        else if (path.endsWith(".lst")) nLST++;
        else if (path.endsWith(".plt")) nPLT++;
        else if (path.endsWith(".b64")) nB64++;
        else if (path.endsWith(".md5")) nMD5++;
        else if (path.endsWith(".mol")) nMOL++;

        if (path.endsWith(".gif")  ||  path.endsWith(".png"))
        {
            if (path.endsWith(".gif")) response.setContentType("image/gif");   // .gif
            else                       response.setContentType("image/png");   // .png
            ServletOutputStream out = response.getOutputStream();
            sendBinaryOutput(out, path, smiles, mol, reactionSmiles, md5, cache,
                             width, height, hasBorder, backColor, transparent != null, labels, title, footer, forColorBlind);
            return;
        } 
        else
        {
            if (path.endsWith(".hex")  || path.endsWith(".bin")  || path.endsWith(".lst")  ||
                path.endsWith(".plt")  || path.endsWith(".b64")  ||
                path.endsWith(".md5"))
                response.setContentType("text/plain");
            else if (path.endsWith(".mol"))
                response.setContentType("chemical/x-mdl-molfile");
            else
                response.setContentType("text/html");   // all other cases
            PrintWriter out = response.getWriter();
            sendTextOutput(out, path, smiles, molfile, mol, reactionSmiles, nbits, asQuery, md5,
                           width, height, hasBorder, backColor, transparent != null, labels, title, footer,
                           request.getRequestURL().toString(), forColorBlind);
        }
    }

    /**
     * Copies the "in" reader to the "out" writer.
     *
     * @param in a BufferedReader
     * @param out a PrintWriter
     */
    private void copyFile(BufferedReader in, PrintWriter out)
            throws IOException {
        int chars;
        char[] buf = new char[4096];
        while ((chars = in.read(buf, 0, buf.length)) != -1) {
            out.write(buf, 0, chars);
            out.flush();
        }
    }

    void sendGif(MOLCanvas canvas, String reactionSmiles,
                 OutputStream output, int width, int height, String[] labels, String title, String footer, boolean forColorBlind)
    {
        // encode the image of the canvas as a GIF
        Graphics2D g2 = null;
        try
        {
            BufferedImage image =
                new BufferedImage(width,height,BufferedImage.TYPE_BYTE_INDEXED);
            g2 = image.createGraphics();
            g2.setColor(Color.white);
            g2.fillRect(0,0,width,height);
            if (reactionSmiles != null)
            {
                String rsmiles = reactionSmiles.substring(0, reactionSmiles.indexOf('>'));
                String psmiles = reactionSmiles.substring(reactionSmiles.lastIndexOf('>') + 1);
                // System.err.println("createDepictionCommands: rsmiles = "+rsmiles);
                // System.err.println("createDepictionCommands: psmiles = "+psmiles);

                StandardMOL rh = null, ph = null;
                synchronized (Mol2ImageServlet.class)
                {
                    String molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(rsmiles, 0);
                    BufferedReader br = new BufferedReader(new StringReader(molfile));
                    FortranInputStream in = new FortranInputStream(br);
                    rh = new StandardMOL();
                    rh.readMOLFile(in);
                    br.close();
                    rh.recolor();
                    rh.clearTopography();
                    molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(psmiles, 0);
                    br = new BufferedReader(new StringReader(molfile));
                    in = new FortranInputStream(br);
                    ph = new StandardMOL();
                    ph.readMOLFile(in);
                    br.close();
                    ph.recolor();
                    ph.clearTopography();
                }

                // Molecule r = SmilesParser.smilesToMolecule(rsmiles);
                // r.recolor();
                // Molecule rh;
                // rh = Layout.layoutMolecule(r);
                // rh.clearTopography();
                // Molecule p = SmilesParser.smilesToMolecule(psmiles);
                // p.recolor();
                // Molecule ph;
                // ph = Layout.layoutMolecule(p);
                // ph.clearTopography();
                Molecule.perceiveReaction(rh, ph);
                String[] reactant_commands = new MoleculeDepicter().computeDepiction(rh, 0, 0,
                                                                                     MoleculeDepicter.USE_COLORS, (int [])null, (String[]) null);
                String[] product_commands = new MoleculeDepicter().computeDepiction(ph, 0, 0,
                                                                                    MoleculeDepicter.USE_COLORS, (int [])null, (String[]) null);
                Dimension rdim = MetaDraw.preferredPlotSize(reactant_commands);
                Dimension pdim = MetaDraw.preferredPlotSize(product_commands);
                Rectangle box = new Rectangle((int) (width * 2 / 5), height);
                MetaDraw.depictCommands(g2, reactant_commands, box, true, forColorBlind);
                box = new Rectangle((int) (width * 3 / 5), 0,
                                    (int) (width * 2 / 5), height);
                MetaDraw.depictCommands(g2, product_commands, box, true, forColorBlind);
                // draw reaction arrow
                g2.drawLine((int) (width * 21 / 50), (int) (height * 25 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
                g2.drawLine((int) (width * 28 / 50), (int) (height * 24 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
                g2.drawLine((int) (width * 28 / 50), (int) (height * 26 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
            }
            else
            {
                String[] commands = canvas.getCommands();
                MetaDraw.depictCommands(g2,
                                        commands,
                                        new Rectangle(width,height),
                                        true,
                                        forColorBlind);
            }
            if (labels != null)
            {
                Font font = g2.getFont();
                FontRenderContext frc = g2.getFontRenderContext();
                int offset = -1;
                for (int i=0; i<labels.length; i++)
                {
                    Rectangle2D rect = font.getStringBounds(labels[i], frc);
                    ((Graphics)g2).drawString(labels[i], 2, height+offset);
                    offset -= (int)rect.getHeight();
                }
            }

            if (title != null &&  title.trim().length() > 0)
            {
                // remenmber old font
                Font font = g2.getFont();
                // set a bigger font for titles
                g2.setFont(font.deriveFont((float)(Math.min(12.0,font.getSize2D()*2.5))));
                FontRenderContext frc = g2.getFontRenderContext();
                // get bounding box
                Rectangle2D rect = g2.getFont().getStringBounds(title, frc);
                ((Graphics)g2).drawString(title, (int)(width/2-rect.getWidth()/2), (int)(3*rect.getHeight()/2));
                // restore old font
                g2.setFont(font);
            }
            if (footer != null &&  footer.trim().length() > 0)
            {
                // remenmber old font
                Font font = g2.getFont();
                // set a bigger font for footers
                g2.setFont(font.deriveFont((float)(Math.min(10.0,font.getSize2D()*1.6))));
                FontRenderContext frc = g2.getFontRenderContext();
                // get bounding box
                Rectangle2D rect = g2.getFont().getStringBounds(footer, frc);
                ((Graphics)g2).drawString(footer, (int)(width/2-rect.getWidth()/2), (int)(height-rect.getHeight()/2));
                // restore old font
                g2.setFont(font);
            }

            g2.dispose(); g2 = null;
            ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
            ImageIO.write(image, "gif", outputStream);
            byte[] bytes = outputStream.toByteArray();
            output.write(bytes);
            // ImageIO.write(image, "gif", output);
            // output.flush();
        }
        catch (IOException e)
        {
            // don't log ClientAbortException instances
            if (0 > e.toString().indexOf("Abort"))
                System.err.println("sendGif(): Caught exception: " + e);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        catch (Throwable t)
        {
            t.printStackTrace();
        }
        finally
        {
            if (g2 != null)
            {
                g2.dispose();
                System.err.println("GIF: Out of bounds dispose()");
            }
        }
    }

    byte[] getPngBytes(MOLCanvas canvas,
                       String reactionSmiles,
                       int width, int height,
                       boolean hasBorder,
                       Color backColor,
                       boolean transparent,
                       String[] labels,
                       String title,
                       String footer,
                       boolean forColorBlind)
    {
        // encode the image of the canvas as a PNG
        byte[] result = new byte[0];
        Graphics2D g2 = null;
        try
        {
            BufferedImage image = null;
            if (transparent)
                image = new BufferedImage(width,height,BufferedImage.TYPE_INT_ARGB);
            else
                image = new BufferedImage(width,height,BufferedImage.TYPE_BYTE_INDEXED);
            g2 = image.createGraphics();
            g2.setColor(backColor);
            g2.fillRect(0,0,width,height);
            if (hasBorder)
            {
                g2.setColor(Color.BLACK);
                g2.drawRect(0, 0, width-1, height-1);
            }
            if (reactionSmiles != null)
            {
                String rsmiles = reactionSmiles.substring(0, reactionSmiles.indexOf('>'));
                String psmiles = reactionSmiles.substring(reactionSmiles.lastIndexOf('>') + 1);
                // System.err.println("createDepictionCommands: rsmiles = "+rsmiles);
                // System.err.println("createDepictionCommands: psmiles = "+psmiles);

                StandardMOL rh = null, ph = null;
                synchronized (Mol2ImageServlet.class)
                {
                    String molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(rsmiles, 0);
                    BufferedReader br = new BufferedReader(new StringReader(molfile));
                    FortranInputStream in = new FortranInputStream(br);
                    rh = new StandardMOL();
                    rh.readMOLFile(in);
                    br.close();
                    rh.recolor();
                    rh.clearTopography();
                    molfile = JNISmi2Mol.getSmi2Mol().smiToMOLWithFlags(psmiles, 0);
                    br = new BufferedReader(new StringReader(molfile));
                    in = new FortranInputStream(br);
                    ph = new StandardMOL();
                    ph.readMOLFile(in);
                    br.close();
                    ph.recolor();
                    ph.clearTopography();
                }
                // Molecule r = SmilesParser.smilesToMolecule(rsmiles);
                // r.recolor();
                // Molecule rh;
                // rh = Layout.layoutMolecule(r);
                // rh.clearTopography();
                // Molecule p = SmilesParser.smilesToMolecule(psmiles);
                // p.recolor();
                // Molecule ph;
                // ph = Layout.layoutMolecule(p);
                // ph.clearTopography();
                Molecule.perceiveReaction(rh, ph);
                String[] reactant_commands = new MoleculeDepicter().computeDepiction(rh, 0, 0,
                                                                                     MoleculeDepicter.USE_COLORS, (int [])null, (String[]) null);
                String[] product_commands = new MoleculeDepicter().computeDepiction(ph, 0, 0,
                                                                                    MoleculeDepicter.USE_COLORS, (int [])null, (String[]) null);
                // System.err.println("reactant commands = " + reactant_commands);
                // System.err.println("product commands = " + product_commands);
                Dimension rdim = MetaDraw.preferredPlotSize(reactant_commands);
                Dimension pdim = MetaDraw.preferredPlotSize(product_commands);
                Rectangle box = new Rectangle((int) (width * 2 / 5), height);
                MetaDraw.depictCommands(g2, reactant_commands, box, true, forColorBlind);
                box = new Rectangle((int) (width * 3 / 5), 0,
                                    (int) (width * 2 / 5), height);
                MetaDraw.depictCommands(g2, product_commands, box, true, forColorBlind);
                // draw reaction arrow
                g2.drawLine((int) (width * 21 / 50), (int) (height * 25 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
                g2.drawLine((int) (width * 28 / 50), (int) (height * 24 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
                g2.drawLine((int) (width * 28 / 50), (int) (height * 26 / 50),
                            (int) (width * 29 / 50), (int) (height * 25 / 50));
            }
            else
            {
                String[] commands = canvas.getCommands();
                MetaDraw.depictCommands(g2,
                                        commands,
                                        new Rectangle(width,height),
                                        true,
                                        forColorBlind);
            }
            if (labels != null)
            {
                Font font = g2.getFont();
                FontRenderContext frc = g2.getFontRenderContext();
                int offset = 0;
                for (int i=0; i<labels.length; i++)
                {
                    Rectangle2D rect = font.getStringBounds(labels[i], frc);
                    offset -= (int)rect.getHeight();
                    ((Graphics)g2).drawString(labels[i], 2, height+offset);
                    offset--;
                }
            }
            if (title != null)
            {
                // remenmber old font
                Font font = g2.getFont();
                // set a bigger font for titles
                g2.setFont(font.deriveFont((float)(Math.min(12.0,font.getSize2D()*2.5))));
                FontRenderContext frc = g2.getFontRenderContext();
                // get bounding box
                Rectangle2D rect = g2.getFont().getStringBounds(title, frc);
                ((Graphics)g2).drawString(title, (int)(width/2-rect.getWidth()/2), (int)(3*rect.getHeight()/2));
                // restore old font
                g2.setFont(font);
            }
            if (footer != null &&  footer.trim().length() > 0)
            {
                // remenmber old font
                Font font = g2.getFont();
                // set a bigger font for footers
                g2.setFont(font.deriveFont((float)(Math.min(10.0,font.getSize2D()*1.6))));
                FontRenderContext frc = g2.getFontRenderContext();
                // get bounding box
                Rectangle2D rect = g2.getFont().getStringBounds(footer, frc);
                ((Graphics)g2).drawString(footer, (int)(width/2-rect.getWidth()/2), (int)(height - rect.getHeight()/2));
                // restore old font
                g2.setFont(font);
            }
            g2.dispose(); g2 = null;
            ByteArrayOutputStream output = new ByteArrayOutputStream();
            if (transparent)
            {
                final BufferedImage transparentImage = imageToBufferedImage(makeColorTransparent(image, backColor));
                ImageIO.write(transparentImage, "png", output);
                // output.flush();
                result = output.toByteArray();
            }
            else
            {
                ImageIO.write(image, "png", output);
                // output.flush();
                result = output.toByteArray();
            }
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        finally
        {
            if (g2 != null)
            {
                g2.dispose();
                System.err.println("PNG: Out of bounds dispose()");
            }
        }
        return result;
    }

   /**
    * Convert Image to BufferedImage.
    *
    * @param image Image to be converted to BufferedImage.
    * @return BufferedImage corresponding to provided Image.
    */
   private static BufferedImage imageToBufferedImage(final Image image)
   {
      final BufferedImage bufferedImage =
         new BufferedImage(image.getWidth(null), image.getHeight(null), BufferedImage.TYPE_INT_ARGB);
      final Graphics2D g2 = bufferedImage.createGraphics();
      g2.drawImage(image, 0, 0, null);
      g2.dispose();
      return bufferedImage;
    }

   /**
    * Make provided image transparent wherever color matches the provided color.
    *
    * @param im BufferedImage whose color will be made transparent.
    * @param color Color in provided image which will be made transparent.
    * @return Image with transparency applied.
    *
    * See: http://www.javaworld.com/community/node/7629
    */
   public static Image makeColorTransparent(final BufferedImage im, final Color color)
   {
      final ImageFilter filter = new RGBImageFilter()
      {
         // the color we are looking for (white)... Alpha bits are set to opaque
         public int markerRGB = color.getRGB() | 0xFFFFFFFF;

         public final int filterRGB(final int x, final int y, final int rgb)
         {
            if ((rgb | 0xFF000000) == markerRGB)
            {
               // Mark the alpha bits as zero - transparent
               return 0x00FFFFFF & rgb;
            }
            else
            {
               // nothing to do
               return rgb;
            }
         }
      };

      final ImageProducer ip = new FilteredImageSource(im.getSource(), filter);
      return Toolkit.getDefaultToolkit().createImage(ip);
   }

    void sendPng(MOLCanvas canvas,
                 String reactionSmiles,
                 OutputStream output,
                 int width, int height,
                 boolean hasBorder,
                 Color backColor,
                 boolean transparent,
                 String[] labels,
                 String title,
                 String footer,
                 boolean forColorBlind)
    {
        // encode the image of the canvas as a PNG
        try
        {
            output.write(getPngBytes(canvas, reactionSmiles, width, height, hasBorder, backColor, transparent, labels, title, footer, forColorBlind));
            // output.flush();
        }
        catch (IOException e)
        {
            // don't log ClientAbortException instances
            if (0 > e.toString().indexOf("Abort"))
                System.err.println("sendPng(): Caught exception: " + e);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    static String bytesToHex(byte[] bytes)
    {
            StringBuffer hex = new StringBuffer(bytes.length * 2);
            int intVal;
            for (int i = 0; i < bytes.length; i++) {
                    intVal = 0xFF & bytes[i];
                    if (intVal < 16)
                            hex.append('0');
                    hex.append(Integer.toHexString(intVal));
            }
            return hex.toString();
    }

    /**
     * Implements a generic LRU cache.
     */
    static class MOLCache<K, V> extends LinkedHashMap<K, V>
    {
        private static int MAX_ENTRIES = 500;
        public MOLCache(int initialCapacity, float loadFactor, boolean accessOrder)
        {
            super(initialCapacity, loadFactor, accessOrder);
            if (MAX_ENTRIES < 5* initialCapacity) MAX_ENTRIES = 5*initialCapacity;
        }

        public synchronized V put(K key, V pngBytes)
        {
            return super.put(key, pngBytes);
        }

        public synchronized V get(String key)
        {
            return super.get(key);
        }


        protected boolean removeEldestEntry(Map.Entry eldest)
        {
           return size() > MAX_ENTRIES;
        }
        
    }

    private void sendBinaryOutput(OutputStream out, String path, String smiles, MOLCanvas mol, String reactionSmiles,
                                  String md5, MOLCache<String, byte[]> cache,
                                  int width, int height, boolean hasBorder, Color backColor, boolean isTransparent, String[] labels,
                                  String title, String footer, boolean forColorBlind)
        throws IOException
    {
        // handle lookup case quickly
        if (path.endsWith(".png")  &&  md5 != null)
        {
            byte[] md5Bytes = cache.get(md5);
            if (md5Bytes != null)
            {
                out.write(md5Bytes);
                return;
            }
        }
        else if (path.endsWith(".gif"))
        {
            sendGif(mol, reactionSmiles, out, width, height, labels, title, footer, forColorBlind);
            return;
        }
        else if (path.endsWith(".png"))
        {
            if (md5 != null)
            {
                byte[] md5Bytes = cache.get(md5);
                if (md5Bytes != null) out.write(md5Bytes);
            }
            else
                sendPng(mol, reactionSmiles, out, width, height, hasBorder, backColor, isTransparent, labels, title, footer, forColorBlind);
            return;
        }
    }

    private void sendTextOutput(PrintWriter out, String path, String smiles, String molfile, MOLCanvas mol, String reactionSmiles,
                                int nbits, boolean asQuery,
                                String md5, int width, int height, boolean hasBorder, Color backColor, boolean isTransparent,
                                String[] labels, String title, String footer,
                                String url, boolean forColorBlind)
        throws IOException
    {
        // generate fingerprint using text formats
        if (path.endsWith(".hex")  ||  path.endsWith(".bin")  ||  path.endsWith(".lst")  ||  path.endsWith(".mol"))
        {
            for (int i=3; i<15; i++)    // make sure FP length is power of 2
                if (nbits < (1<<i))
                {
                    nbits = (1<<i);
                    break;
                }
            byte[] molfp = null;
            synchronized (Mol2ImageServlet.class)
            {
                JNIMatch matcher = JNIMatch.getMatcher();
                if (molfile == null  &&  smiles != null  && !smiles.equals(""))
                molfile = JNISmi2Mol.getSmi2Mol().smiToMOL(smiles);
                // byte[] molfp = new byte[nbits/8];
                molfp = matcher.getFingerprintFromCT(molfile, asQuery, JNIMatch.USE_ALL_BITS, nbits/8);
            }
            if (path.endsWith(".hex"))
                out.println(bytesToHex(molfp));
            else if (path.endsWith(".mol"))
            {
                if (molfile == null || molfile.trim().equals(""))
                {
                    molfile = // No Structure MOL File
                      "\n" +
                      "  -NONE-  12070108312D\n" +
                      "\n" +
                      "  0  0  0  0  0  0  0  0  0  0999 V2000\n" +
                      "M  END\n";
                }
                out.println(molfile);
            }
            else if (path.endsWith(".bin"))
            {
                for (int i=0; i<molfp.length; i++)
                    for (int j=0; j<8; j++)
                        if ((molfp[i]&(byte)(1<<j)) != 0) out.print('1');
                        else                              out.print('0');
                out.println();
            }
            else // if (path.endsWith(".lst"))
            {
                boolean first = true;
                for (int i=0; i<molfp.length; i++)
                    for (int j=0; j<8; j++)
                        if ((molfp[i]&(byte)(1<<j)) != 0)
                        {
                            if (first) out.print(8*i+j);
                            else out.print(","+(8*i+j));
                            first=false;
                        }
                out.println();
            }
            return;
        }
        else if (path.endsWith(".plt"))
        {
            String[] commands = mol.getCommands();
            for (int i=0; i<commands.length; i++)
                out.println(commands[i]);
            return;
        }
        else if (path.endsWith(".b64"))        // return base64 coded picture
        {
            byte[] pngBytes = getPngBytes(mol, reactionSmiles,
                                          width, height, hasBorder, backColor, isTransparent, labels, title, footer, forColorBlind);
            out.print(Base64.encodeBytes(pngBytes));
        }
        else if (path.endsWith(".md5")  ||     // cache a PNG and return MD5 hash
                 path.endsWith(".html"))       // use cached PNG to populate an image tag in the HTML output
        {
            byte[] pngBytes = getPngBytes(mol, reactionSmiles,
                                          width, height, hasBorder, backColor, isTransparent, labels, title, footer, forColorBlind);
            try
            {
                MessageDigest md = MessageDigest.getInstance("MD5");
                md.reset();
                md.update(pngBytes);
                byte[] digest = md.digest();
                StringBuffer hexString = new StringBuffer();
                for (int i=0;i<digest.length;i++)
                    hexString.append(Integer.toHexString(0xFF & digest[i]));
                md5 = hexString.toString();
                cache.put(md5, pngBytes);
                if (path.endsWith(".md5"))
                {
                    out.println(md5);
                }
                else
                {
                    out.println("<html><body>");
                    out.println("<img src='" + url.replaceFirst(".html",".png") + "?md5=" + md5 + "'>");
                    out.println("</body></html>");
                }
            }
            catch (Exception e)
            {
                e.printStackTrace();
            }
            return;
        }
        else
        {
            if (path.endsWith("/")  ||  path.contains("version"))
            {
                out.print("version number = 1.2.0");
                out.println(", build number = " + buildNumber);
            }
            else
            {
                out.print("unknown path");
            }
        }
    }

    // Uses simple String manipulations on SMILES (no JNI calls)
    static String abbreviateChains(String smiles)
    {
        if (smiles == null  || smiles.length() < 5) return smiles;
        boolean changed;
        // 'after BOS before (' case
        do
        {
            changed = false;

            Pattern p = Pattern.compile("^(C{5,99})C\\(");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+1;
                int nC = group.length();
                smiles = m.replaceFirst("{C"+nC+"H"+nH+"}C(");
                changed = true;
            }
        } while (changed);
        // 'after . before (' case
        do
        {
            changed = false;

            Pattern p = Pattern.compile("\\.(C{5,99})C\\(");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+1;
                int nC = group.length();
                smiles = m.replaceFirst(".{C"+nC+"H"+nH+"}C(");
                changed = true;
            }
        } while (changed);
        // 'after BOS before non-C' case
        do
        {
            changed = false;

            Pattern p = Pattern.compile("^(C{5,99})([NOS])");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+1;
                int nC = group.length();
                smiles = m.replaceFirst("{C"+nC+"H"+nH+"}"+m.group(2));
                changed = true;
            }
        } while (changed);
        // 'after . before non-C' case
        do
        {
            changed = false;

            Pattern p = Pattern.compile("\\.(C{5,99})([NOS])");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+1;
                int nC = group.length();
                smiles = m.replaceFirst(".{C"+nC+"H"+nH+"}"+m.group(2));
                changed = true;
            }
        } while (changed);
        // end of parenthesis case
        do
        {
            changed = false;

            Pattern p = Pattern.compile("(C{5,99})C\\)");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+3;
                int nC = group.length()+1;
                smiles = m.replaceFirst("{C"+nC+"H"+nH+"})");
                changed = true;
            }
        } while (changed);
        // end of string case
        do
        {
            changed = false;
            Pattern p = Pattern.compile("(C{5,99})C$");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+3;
                int nC = group.length()+1;
                smiles = m.replaceFirst("{C"+nC+"H"+nH+"}");
                changed = true;
            }
        } while (changed);
        // end of string before '.' case
        do
        {
            changed = false;
            Pattern p = Pattern.compile("(C{5,99})C\\.");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(1);
                int nH = 2*group.length()+3;
                int nC = group.length()+1;
                smiles = m.replaceFirst("{C"+nC+"H"+nH+"}.");
                changed = true;
            }
        } while (changed);
        // PEG case
        do
        {
            changed = false;
            Pattern p = Pattern.compile("(CCO){4,99}C");
            Matcher m = p.matcher(smiles);
            if (m.find())
            {
                String group = m.group(0);
                int n = (group.length()-1)/3;
                smiles = m.replaceFirst("{PEG"+n+"}");
                changed = true;
            }
        } while (changed);
        return smiles;
    }

    /**
     * Use this main() method to run as stand-alone.
     */
    public static void main(String[] argv)
    {
    }
}
