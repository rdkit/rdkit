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
This directory contains the following executables and library modules:


canonizer.exe
-------------
Canonicalizes input SMILES (with or w/o IDs). Run

    canonizer -h

to get a (short) description of the command line parameters.


mol2smi.exe
-----------
Converts an input SDF or MOL file to SMILES. Run

    mol2smi -h

to get a (short) description of the command line parameters.


mol2tbl.exe
-----------
A variant of mol2smi that works on SDF files and preserves selected data fields
in the output. Run

    mol2tbl -h

to get a (short) description of the command line parameters.


smi2mol.exe
-----------
Computes a MOL file from a SMILES by inventing atom coordinates. Run

    smi2mol -h

to get a (short) description of the command line parameters.


tomcat_lib/*
------------
The contents of this directory needs to be put on the BOOTCLASSPATH of a Tomcat
server to provide the JNI binaries and the glue code for a service that would
generate pictures from MOL files and SMILES (among a few other things like
fingerprints). The reason for putting these files on the BOOTCLASSPATH is that
they should only be loaded once.


depicter.war
------------
This is the Web application that provided the depiction service when deployed
on a Tomcat server that has the above JNI bits installed.

The default file index.html describes (some) of the URL parameters and the JSP
SmilesVisualizer.jsp is an interactive renderer of SMILES typed into an input
field.

struchk.exe
-----------
Efficient augmented atom-based chemical structure standardization tool. It uses
the files checkfgs.* as input. Run

    struchk -h

to get a (short) description of the command line parameters.

depict32.dll
------------
Dynamic link library used to render pictures from SMILES and MOL files.


