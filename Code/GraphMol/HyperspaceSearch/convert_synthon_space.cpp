//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is a simple program that takes a synthon file in Chemspace/Enamine
// format and converts it into the binary format used for searching.  It
// changes the connector labels from [U], [Np], [Pu], [Am] to [1*], [2*],
// [3*], [4*] along the way.
// The format of the input file can be either (tab-separated, release column)
// optional.
// SMILES	synton_id	synton#	reaction_id release
// [1*]Nc1c([2*])cccc1	1-1	0	doebner-miller-quinoline  3
// [1*]Nc1cc(C)ccc1[2*]	1-2	0	doebner-miller-quinoline  3
// [1*]Nc1ccc(C(Cl)(Cl)Cl)cc1[2*]	1-3	0	doebner-miller-quinoline
// 3
//  etc.
//
// or
// SMILES,synton_id,synton_role,reaction_id
// O=C(c1ccccc1C(F)(F)F)N1CCCN([U])CC1,236434,synton_1,a2
// [U]c1ncc2c(n1)CCCC2,250139,synton_2,a2
// CNC(=O)c1ccnc([U])c1,19571,synton_2,a2
// etc.
//
// It takes 2 arguments, the names of the input and output files.

#include <iostream>

#include <GraphMol/HyperspaceSearch/Hyperspace.h>

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "convertSynthonSpace takes 2 arguments, the name of the text"
              << " file and the name of the binary output file." << std::endl;
    return (1);
  }
  std::string inFile(argv[1]);
  std::string outFile(argv[2]);
  std::cout << "Converting " << inFile << " to " << outFile << std::endl;
  RDKit::HyperspaceSearch::Hyperspace hyperspace;
  hyperspace.readTextFile(inFile);
  hyperspace.summarise(std::cout);
  hyperspace.writeDBFile(outFile);
}