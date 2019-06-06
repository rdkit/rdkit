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
#include <GraphMol/RDKitBase.h>

class JSMol {
 public:
  JSMol() : d_mol(nullptr){};
  JSMol(RDKit::ROMol *mol) : d_mol(mol){};
  std::string get_smiles() const;
  std::string get_molblock() const;
  std::string get_inchi() const;
  std::string get_svg() const;
  std::string get_svg_with_highlights(const std::string &details) const;
  std::string get_substruct_match(const JSMol &q) const;
  std::string get_descriptors() const;
  std::string get_morgan_fp(unsigned int radius, unsigned int len) const;
  std::string get_morgan_fp() const { return get_morgan_fp(2, 2048); };
  bool is_valid() const { return d_mol.get() != nullptr; };

 private:
  std::unique_ptr<RDKit::ROMol> d_mol;
};

std::string get_inchikey_for_inchi(const std::string &input);
JSMol *get_mol(const std::string &input);
JSMol *get_qmol(const std::string &input);
std::string version();
