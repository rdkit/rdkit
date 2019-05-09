//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>

std::string get_smiles(const std::string &input);
std::string get_inchi(const std::string &input);
std::string get_inchikey_for_inchi(const std::string &input);
std::string get_svg(const std::string &input);
std::string get_pkl(const std::string &input);
std::string version();
int ping();
