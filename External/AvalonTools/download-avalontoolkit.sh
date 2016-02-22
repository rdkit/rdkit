#!/bin/bash
#
#  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# This script will download InChI software distribution version 1.04 and 
# places it (temporarily) in a /tmp directory.  This directory will 
# automatically dissappear after you reboot the machine.

if ! which wget > /dev/null
then
	echo "**Error: wget not found!"
	exit 1
fi

DIR=`pwd`
TEMPDIR=`mktemp -d -t rdkit-avalontoolkit-XXX`
if [[ $DIR =~ External/AvalonTools$ ]]
then
	mkdir -p src
	echo "================================================================"
	echo "Downloading AvalonTools software distribution version 1.04"
	echo "  http://sourceforge.net/projects/avalontoolkit/files/AvalonToolkit_1.2/AvalonToolkit_1.2.0.source.tar"
	echo "  ====>"
	echo "  $TEMPDIR"
	echo "================================================================"
	cd $TEMPDIR
	wget http://sourceforge.net/projects/avalontoolkit/files/AvalonToolkit_1.2/AvalonToolkit_1.2.0.source.tar
	
	echo "================================================================"
	echo "Unarchiving"
	echo "================================================================"
	tar xvf AvalonToolkit_1.2.0.source.tar
	echo "================================================================"
	echo "Applying AvalonToolkit_1.2_patch.txt"
	echo "================================================================"
  patch -p0 < $DIR/AvalonToolkit_1.2_patch.txt
	echo "================================================================"
	echo "Copying files"
	echo "================================================================"
	cp -R SourceDistribution "$DIR"
	echo "================================================================"
	echo "Removing temporary files"
	echo "================================================================"
	cd $DIR
	rm -rf $TEMPDIR
	echo "================================================================"
	echo "Done!"
	echo "Make sure you (re)run cmake before running make"
	echo "================================================================"
else
	echo '**Error: you must invoke this script from within the directory'
	echo '         $RDKIT_SOURCE_ROOT/External/AvalonTools'
	exit 1
fi
