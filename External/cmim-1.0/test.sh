#!/bin/sh

##############################################################################
# This program is free software; you can redistribute it and/or              #
# modify it under the terms of the GNU General Public License                #
# version 2 as published by the Free Software Foundation.                    #
#                                                                            #
# This program is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of                 #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          #
# General Public License for more details.                                   #
#                                                                            #
# Written by François Fleuret                                                #
# Contact <francois.fleuret@epfl.ch> for comments & bug reports              #
# Copyriche (C) 2004 EPFL                                                    #
##############################################################################

# $Id: test.sh,v 1.1 2005/03/03 15:52:35 fleuret Exp $

./create_samples

for fs in cmim mim random; do
    ./cmim --nb-features 10 --feature-selection $fs --classifier bayesian --train ./train.dat /tmp/$fs.clf
    ./cmim --silent --test /tmp/$fs.clf test.dat results.dat | grep ERROR
done
