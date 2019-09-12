# This file is part of the RingDecomposerLib, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2016
# University of Hamburg, ZBH - Center for Bioinformatics
# Niek Andresen, Florian Flachsenberg, Matthias Rarey
# 
# Please cite:
# 
# Kolodzik, A.; Urbaczek, S.; Rarey, M.
# Unique Ring Families: A Chemically Meaningful Description
# of Molecular Ring Topologies.
# J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
# 
# Flachsenberg, F.; Andresen, N.; Rarey, M.
# RingDecomposerLib: An Open-Source Implementation of
# Unique Ring Families and Other Cycle Bases.
# J. Chem. Inf. Model., 2017, 57 (2), pp 122-126

## @file simple_demo.py
## @brief This file containes simple demo for using
## the Python wrapper of the RingDecomposerLib.
##
## usage:
##
##     simple_demo.py
##
## The program takes no parameters.
## It calculates the ring topology for an example graph
## and prints it to stdout.

from __future__ import print_function
import sys
import py_rdl


def calculate_and_print_urfs(edges):
    print('Extensive example')
    print('*' * 80)
    a = py_rdl.Graph.from_edges(edges)

    a.add_edge(edges[-1])

    d = py_rdl.Calculator(a)

    try:
        d.get_nof_urf()
    except py_rdl.RDLError:
        print('caught exception, because calculate wasn\'t called')

    d.calculate()
    print(d.get_nof_urf())

    for urf in d:
        print(urf)
        print(urf.nodes)
        print(urf.edges)
        print(d.get_relevant_cycles_for_urf(urf))

    print('MCB')
    print(d.get_sssr())

    print('RCs')
    print(d.get_relevant_cycles())

    print('RCPs')
    print(d.get_relevant_cycle_prototypes())

def minimal_example(edges):
    print('Minimal example')
    print('*' * 80)
    urf_data = py_rdl.Calculator.get_calculated_result(edges)
    for urf in urf_data:
        print(urf)

if __name__ == '__main__':

    edges = [
        ("0", "1"),
        ("0", "2"),
        ("1", "3"),
        ("2", "3"),
        ("3", "4"),
        ("4", "5"),
        ("4", "6"),
        ("5", "7"),
        ("6", "7"),
        ("7", "0")
    ]

    minimal_example(edges)
    print()
    calculate_and_print_urfs(edges)
