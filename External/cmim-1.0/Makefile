
##############################################################################
# This program is free software; you can redistribute it and/or		     #
# modify it under the terms of the GNU General Public License		     #
# version 2 as published by the Free Software Foundation.		     #
# 									     #
# This program is distributed in the hope that it will be useful, but	     #
# WITHOUT ANY WARRANTY; without even the implied warranty of		     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU	     #
# General Public License for more details.				     #
# 									     #
# Written by François Fleuret                                                #
# Contact <francois.fleuret@epfl.ch> for comments & bug reports              #
# Copyright (C) 2004 EPFL                                                    #
##############################################################################

ifeq ($(CXX),)
  CXX = g++-2.95
endif

ifeq ($(DEBUG),yes)
 LDFLAGS = -lm
 CXXFLAGS = -Wall -g -Wall -DDEBUG -DVERBOSE
else
 LDFLAGS = -lm
 CXXFLAGS =  -ffast-math -fcaller-saves -finline-functions -funroll-all-loops -Wall -O3
 CXXFLAGS =  -Wall
endif

all: cmim.exe create_samples Makefile.depend

test: all
	./test.sh

cmim.exe: misc.o fastentropy.o classifier.o cmim.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

libfastentropy.a: fastentropy.o misc.o
	ar -r $@ $^

create_samples: create_samples.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	\rm -f *.o cmim.exe

Makefile.depend: *.h *.cc Makefile
	$(CC) -M *.cc > Makefile.depend

ARCHIVE_NAME =  cmim-1.0
DIRNAME = cmim-1.0

archive:
	cd ..; tar zcvf ${ARCHIVE_NAME}.tgz \
	                ${DIRNAME}/*.{cc,h,sh,txt} ${DIRNAME}/Makefile ${DIRNAME}/README

-include Makefile.depend
