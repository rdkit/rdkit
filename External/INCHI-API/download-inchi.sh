#!/bin/bash

# This script will download InChI software distribution version 1.03 and place
# in the current directory

if ! which wget > /dev/null
then
	echo "**Error: wget not found!"
	exit 1
fi

DIR=`pwd`
TEMPDIR=`mktemp -d -t rdkit-inchi`
if [[ $DIR =~ External/INCHI-API$ ]]
then
	mkdir -p src
	echo "================================================================"
	echo "Downloading InChI software distribution version 1.03"
	echo "  http://www.iupac.org/inchi/download/version1.03/INCHI-1-API.zip"
	echo "  ====>"
	echo "  $TEMPDIR"
	echo "================================================================"
	cd $TEMPDIR
	wget http://www.iupac.org/inchi/download/version1.03/INCHI-1-API.zip
	echo "================================================================"
	echo "Unarchiving"
	echo "================================================================"
	unzip INCHI-1-API.zip
	echo "================================================================"
	echo "Copying files"
	echo "================================================================"
	cp INCHI-1-API/INCHI_API/inchi_dll/* "$DIR/src"
	echo "================================================================"
	echo "Removing temporary files"
	echo "================================================================"
	cd $DIR
	rm -rf $TEMPDIR
	echo "================================================================"
	echo "Done!"
	echo "================================================================"
else
	echo '**Error: you must invoke this script from within the directory'
	echo '         $RDKIT_SOURCE_ROOT/External/INCHI-API'
	exit 1
fi
