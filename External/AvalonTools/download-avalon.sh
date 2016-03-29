#!/bin/bash
#
# This script will download the currently supported version of the Avalon toolkit
# source and places it (temporarily) in a /tmp directory.  This directory will
# automatically dissappear after you reboot the machine.

if ! which wget > /dev/null
then
	echo "**Error: wget not found!"
	exit 1
fi

FNAME="AvalonToolkit_1.2.0.source.tar"
DIR=`pwd`
TEMPDIR=`mktemp -d -t rdkit-avalon-XXX`
if [[ $DIR =~ External/AvalonTools$ ]]
then
	mkdir -p distrib
	echo "================================================================"
	echo "Downloading Avalon toolkit distribution version 1.2.0"
	echo "  $FNAME"
	echo "  ====>"
	echo "  $TEMPDIR"
	echo "================================================================"
	cd $TEMPDIR
	wget https://sourceforge.net/projects/avalontoolkit/files/AvalonToolkit_1.2/$FNAME

	echo "================================================================"
	echo "Unarchiving"
	echo "================================================================"
	tar xf $FNAME

	echo "================================================================"
	echo "Copying files"
	echo "================================================================"
	cp -r SourceDistribution/* "$DIR/distrib"
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
