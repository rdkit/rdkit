#-------------------------------------------------------------------------
#
# Makefile--
#    Makefile for tutorial
#
# IDENTIFICATION
#    $Header: /cvsroot/pgsql-server/src/tutorial/Makefile,v 1.16 2002/09/05 18:28:46 petere Exp $
#
#-------------------------------------------------------------------------

include $(RDBASE)/Code/rdvars.make

CXXFLAGS = $(BASECXXFLAGS) $(RDKITINCFLAGS) -fpic -I/usr/src/redhat/SOURCES/postgresql-8.0.1/src/include
DLSUFFIX=.so

#DLOBJS= $(RDBASE)/bin/funcs$(DLSUFFIX)
DLOBJS= funcs$(DLSUFFIX)

QUERIES= test1.sql

all: $(DLOBJS) $(QUERIES)

#$(RDBASE)/bin/funcs$(DLSUFFIX): funcs.o
#	g++ -shared -o $(RDBASE)/bin/funcs$(DLSUFFIX) funcs.o
#	g++ -shared -o $(RDBASE)/bin/funcs$(DLSUFFIX) funcs.o $(RDSMILES) $(RDSUBSTRUCT) $(RDKIT) $(RDGENERAL)
funcs$(DLSUFFIX): funcs.o
	g++ -shared -o funcs$(DLSUFFIX) funcs.o  $(RDSMILES) $(RDSUBSTRUCT) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY) $(RDFINGERPRINTS) $(RDSUBGRAPHS) $(VFLIB) $(LAPACKLIB) 

clean:
	rm -f $(DLOBJS)
