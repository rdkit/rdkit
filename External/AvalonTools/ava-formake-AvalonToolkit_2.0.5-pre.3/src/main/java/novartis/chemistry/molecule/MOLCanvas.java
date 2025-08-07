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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import novartis.utilities.FortranInputStream;
import novartis.utilities.MetaDraw;

public class MOLCanvas extends JPanel
{
    String molfile;
    String[] commands;
    public synchronized String[] getCommands()
    {
        return commands;
    }

    boolean new_layout = false;

    boolean colored = true;
    public void setColored(boolean colored) { this.colored = colored; }
    public boolean getColored() { return colored; }

    class MyJLabel extends JLabel
    {
        public MyJLabel(String text)
        {
            super(text);
        }

        public void paintComponent(Graphics g)
        {
            super.paintComponent(g);
        }
    }

    Dimension minBox, myBox;
    MyJLabel noStructure = new MyJLabel("No Structure");

    public MOLCanvas(String molfile, int width, int height, boolean new_layout)
    {
        setOpaque(false);
        myBox = new Dimension(width, height);
        minBox = myBox;

        this.molfile = molfile;
        this.new_layout = new_layout;
        this.setLayout(new BorderLayout());
        noStructure.setHorizontalAlignment(SwingConstants.CENTER);
        noStructure.setVerticalAlignment(SwingConstants.CENTER);
        noStructure.setForeground(Color.black);
        if (molfile == null){
        	this.add(noStructure, BorderLayout.CENTER);
        	noStructure.setVisible(true);
            return; 
        } else
            noStructure.setVisible(false);

        createDepictionCommands(molfile, null, "MOLCanvas");
        invalidate();
    }

    public synchronized void setSize(int width, int height)
    {
        myBox = new Dimension(Math.max(width, minBox.width),
                              Math.max(height, minBox.height));
        repaint();
    }

    private void createDepictionCommands(String molfile,
                                         int[] colors,
                                         String caller)
    {
        try
        {
            BufferedReader br =
                    new BufferedReader(new StringReader(molfile));
            FortranInputStream in = new FortranInputStream(br);

            StandardMOL m = new StandardMOL();
            m.readMOLFile(in);
            br.close();
            m.recolor();
            Molecule mh;
            if (new_layout) mh = Layout.layoutMolecule(m);
            else            mh = m;
            commands =
                new MoleculeDepicter().computeDepiction(mh,
                    0, 0,
                    MoleculeDepicter.USE_COLORS,
                    colors,
                    (String[])null);
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.err.println(caller + ": Caught " + e);
        }
    }

    public synchronized String getMOLFile()
    {
        return molfile;
    }

    public synchronized void setMOLFile(String molfile)
    {
        setMOLFile(molfile, null);
    }

    public synchronized void setMOLFile(String molfile, int[] colors)
    {
        this.molfile = molfile;
        if (molfile != null && !molfile.equals("")) 
        {
            if (molfile.startsWith("sRNA:")) {
            	commands = null;
            	noStructure.setText(molfile.substring(4, molfile.length()));
            	noStructure.setVisible(true);
            } else {
            	createDepictionCommands(molfile, colors, "setMolFile");
            }
            noStructure.setVisible(false);
            invalidate();
        }
         else
        {
           commands = null; // empty canvas
           noStructure.setVisible(true);
           invalidate();
        }
        repaint();
    }

    /**
     * Set the MOLCanvas to a blank JPanel (no 'No Structure').
     */
    public synchronized void setBlank()
    {
       this.molfile = null;
       commands = null; // empty canvas
       noStructure.setVisible(false);
       invalidate();
       repaint();
    }

    public synchronized Dimension getMinimumSize()
    {
        return (minBox);
    }

    public synchronized Dimension getPreferredSize()
    {
        return (minBox);
    }

    public void paint(Graphics g)
    {
        paintComponent(g);
    }

    public synchronized void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        if (commands != null)
        {
            MetaDraw.depictCommands(g,
                                    commands,
                                    new Rectangle(this.getSize()),
                                    colored, false);
        }
        else
        {
            // noStructure.paintComponent(g);
        }
        super.paintBorder(g);
    }

    /**
     * Test driver for MOLCanvas.
     *
     * Creates a JFrame that renders the molecule give on System.in.
     */
    public static void main(String argv[]) throws IOException
    {
        StringBuffer molfile = new StringBuffer();
        if (argv.length == 0) 
        {
            BufferedReader reader =
                new BufferedReader(new InputStreamReader(System.in));
            String line;
            while (null != (line = reader.readLine()))
            {
                molfile.append(line); molfile.append("\n");
            }
        }
        JFrame frame = new JFrame();
        frame.setSize(400, 400);
        MOLCanvas mol = null;
        if (molfile.length() > 0)
            mol = new MOLCanvas(molfile.toString(), 300, 300, false);
        else
            mol = new MOLCanvas(null, 300, 300, false);
        frame.getContentPane().add(mol);
        frame.addWindowListener(
            new WindowAdapter()
            {
                public void windowClosing(WindowEvent e)
                {
                    Window w = e.getWindow();
                    w.setVisible(false);
                    w.dispose();
                    System.exit(0);
                }
            });
        frame.setVisible(true);
    }

}
