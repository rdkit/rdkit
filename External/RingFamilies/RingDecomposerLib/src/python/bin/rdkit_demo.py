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

## @file rdkit_demo.py
## @brief This file containes simple demo for using the
## Python wrapper
## of the RingDecomposerLib together with RDKit.
##
## usage:
##
##     rdkit_demo.py --input <input_file> [--output]
##
## The program takes as input a molecule file
## (.smi or .sdf) and reads in every molecule and
## calculates the ring topology and prints it to stdout.
##
## The parameter `--output` is optional. If specified,
## the script will write a number of SVG images,
## where RCs are marked according to family membership.

from __future__ import absolute_import, print_function

import argparse
import os
import random
import sys

import py_rdl

import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Draw


def analyze_molecule(molecule):
    data = py_rdl.Calculator.get_calculated_result(
               molecule.GetBonds(),              # pass edges
               rdkit.Chem.Bond.GetBeginAtom,     # get first atom
               rdkit.Chem.Bond.GetEndAtom,       # get second atom
               rdkit.Chem.Atom.GetIdx,           # identify atom
               rdkit.Chem.Bond.GetIdx            # identify edge
           )

    print('URFs')
    for urf in data:
        print(urf)
        print(urf.weight)
        print(urf.nodes)
        print(urf.edges)
        print(data.get_relevant_cycles_for_urf(urf))

    print('URF for each atom')
    for a in molecule.GetAtoms():
        print(a.GetIdx(), data.get_urfs_for_node(a.GetIdx()))

    print('URF for each bond')
    for b in molecule.GetBonds():
       print(b.GetIdx(), data.get_urfs_for_edge(b.GetIdx()))

    print('MCB')
    print(data.get_sssr())

    return data


def get_rgb_from_hsv(h, s, v, a):
    def clamp(v, minv, maxv):
        return max(minv, min(v, maxv))

    h = clamp(h, 0.0, 360.0)
    s = clamp(s, 0.0, 1.0)
    v = clamp(v, 0.0, 1.0)

    i = int(h / 60.0)
    r = (h % 60.0) / 60.0
    p = v * (1.0 - s)
    q = v * (1.0 - s * r)
    t = v * (1.0 - s * (1.0 - r))

    if i == 0 or i == 6:
        return (v, t, p, a)
    elif i == 1:
        return (q, v, p, a)
    elif i == 2:
        return (p, v, t, a)
    elif i == 3:
        return (p, q, v, a)
    elif i == 4:
        return (t, p, v, a)
    elif i == 5:
        return (v, p, q, a)

    raise ValueError("Internal color error")


def get_distinct_colors(nof_colors, h=0.0, s=0.8, v=0.4, a=1.0, variation=0.2):
    colors = []

    for i in range(nof_colors):
        h_ = (h + i * 360.0 / nof_colors) % 360.0
        s_ = s + variation * random.random()
        v_ = v + variation + random.random()
        colors.append(get_rgb_from_hsv(h_, s_, v_, a)[:-1])

    return colors


def draw_molecule(mol, filename, higlights):
    rdkit.Chem.AllChem.Compute2DCoords(mol)

    drawer = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(200, 200)
    drawer.DrawMolecule(mol, highlightAtoms=None, highlightBonds=highlights,
        highlightAtomColors=None, highlightBondColors=highlights)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')

    with open(filename, 'w') as outfile:
        outfile.write(svg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('RDKit ring family example')
    parser.add_argument('--input', type=str, required=True,
                        help='Input file (.smi or .sdf)')
    parser.add_argument('--output', action='store_true',
                        help='Write .svg output files with URFs')

    args = parser.parse_args()

    ext = os.path.splitext(args.input)[1]

    if ext == '.smi':
        supplier = rdkit.Chem.SmilesMolSupplier(args.input, titleLine=False)
    elif ext == '.sdf':
        supplier = rdkit.Chem.SDMolSupplier(args.input)
    else:
        raise ValueError('invalid extension: %s' % ext)

    for i, mol in enumerate(supplier):
        print("Molecule", i)
        if mol is None:
            print("ERROR")
            continue
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = "NO_NAME"
        print(name)
        data = analyze_molecule(mol)
        if args.output:
            rcps = data.get_relevant_cycle_prototypes()
            rcf_colors = get_distinct_colors(len(rcps))

            urf_colors = get_distinct_colors(data.get_nof_urf())
            for j, r in enumerate(rcps):
                outfile_name = args.input + '_{}_rcp_{}.svg'.format(i, j)
                highlights = {b: rcf_colors[j] for b in r.edges}
                draw_molecule(mol, outfile_name, highlights)
            for j, r in enumerate(data.get_relevant_cycles()):
                outfile_name = args.input + '_{}_rc_{}_rcf.svg'.format(i, j)
                highlights = {b: rcf_colors[r.rcf.index] for b in r.edges}
                draw_molecule(mol, outfile_name, highlights)

                outfile_name = args.input + '_{}_rc_{}_urf.svg'.format(i, j)
                highlights = {b: urf_colors[r.urf.index] for b in r.edges}
                draw_molecule(mol, outfile_name, highlights)
