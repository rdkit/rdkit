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
  JSMol(RDKit::RWMol *mol) : d_mol(mol){};
  std::string get_smiles() const;
  std::string get_cxsmiles() const;
  std::string get_molblock() const;
  std::string get_v3Kmolblock() const;
  std::string get_inchi() const;
  std::string get_svg(unsigned int width, unsigned int height) const;
  std::string get_svg() const {
    return get_svg(d_defaultWidth, d_defaultHeight);
  };
  std::string get_svg_with_highlights(const std::string &details) const;
  std::string get_substruct_match(const JSMol &q) const;
  std::string get_substruct_matches(const JSMol &q) const;
  std::string get_descriptors() const;
  std::string get_morgan_fp(unsigned int radius, unsigned int len) const;
  std::string get_morgan_fp() const { return get_morgan_fp(2, 2048); };
  std::string condense_abbreviations(double maxCoverage, bool useLinkers);
  std::string condense_abbreviations() {
    return condense_abbreviations(0.4, false);
  };
  std::string condense_abbreviations_from_defs(const std::string &definitions,
                                               double maxCoverage,
                                               bool areLinkers);
  std::string generate_aligned_coords(const JSMol &templateMol,bool useCoordGen);
  std::string generate_aligned_coords(const JSMol &templateMol) { return generate_aligned_coords(templateMol,false);};
  

  bool is_valid() const { return d_mol.get() != nullptr; };

  // functionality primarily useful in ketcher
  std::string get_stereo_tags() const;
  std::string get_aromatic_form() const;
  std::string get_kekule_form() const;
  std::string get_new_coords(bool useCoordGen) const;
  std::string get_new_coords() const { return get_new_coords(false); };
  std::string remove_hs() const;
  std::string add_hs() const;

  std::unique_ptr<RDKit::RWMol> d_mol;
  static constexpr unsigned int d_defaultWidth = 250;
  static constexpr unsigned int d_defaultHeight = 200;
};

std::string get_inchikey_for_inchi(const std::string &input);
JSMol *get_mol(const std::string &input);
JSMol *get_qmol(const std::string &input);
std::string version();
void prefer_coordgen(bool prefer);
