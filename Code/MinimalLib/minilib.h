//
//
//  Copyright (C) 2019-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>

class JSMol {
 public:
  JSMol() : d_mol(nullptr) {}
  JSMol(RDKit::RWMol *mol) : d_mol(mol) {}
  std::string get_smiles() const;
  std::string get_cxsmiles() const;
  std::string get_molblock() const;
  std::string get_v3Kmolblock() const;
  std::string get_inchi() const;
  std::string get_json() const;
  std::string get_svg(unsigned int width, unsigned int height) const;
  std::string get_svg() const {
    return get_svg(d_defaultWidth, d_defaultHeight);
  }
  std::string get_svg_with_highlights(const std::string &details) const;
  std::string get_substruct_match(const JSMol &q) const;
  std::string get_substruct_matches(const JSMol &q) const;
  std::string get_descriptors() const;
  std::string get_morgan_fp(unsigned int radius, unsigned int len) const;
  std::string get_morgan_fp() const { return get_morgan_fp(2, 2048); }
  std::string get_morgan_fp_as_binary_text(unsigned int radius,
                                           unsigned int len) const;
  std::string get_morgan_fp_as_binary_text() const {
    return get_morgan_fp(2, 2048);
  }
  std::string condense_abbreviations(double maxCoverage, bool useLinkers);
  std::string condense_abbreviations() {
    return condense_abbreviations(0.4, false);
  }
  std::string condense_abbreviations_from_defs(const std::string &definitions,
                                               double maxCoverage,
                                               bool areLinkers);
  std::string generate_aligned_coords(const JSMol &templateMol,
                                      bool useCoordGen,
                                      bool allowOptionalAttachments,
                                      bool acceptFailure);
  std::string generate_aligned_coords(const JSMol &templateMol,
                                      bool useCoordGen,
                                      bool allowOptionalAttachments) {
    return generate_aligned_coords(templateMol, useCoordGen,
                                   allowOptionalAttachments, true);
  };
  std::string generate_aligned_coords(const JSMol &templateMol,
                                      bool useCoordGen) {
    return generate_aligned_coords(templateMol, useCoordGen, false, true);
  }
  std::string generate_aligned_coords(const JSMol &templateMol) {
    return generate_aligned_coords(templateMol, false, false, true);
  }

  bool is_valid() const { return d_mol.get() != nullptr; }

  // functionality primarily useful in ketcher
  std::string get_stereo_tags() const;
  std::string get_aromatic_form() const;
  std::string get_kekule_form() const;
  std::string get_new_coords(bool useCoordGen) const;
  std::string get_new_coords() const { return get_new_coords(false); }
  std::string remove_hs() const;
  std::string add_hs() const;

  std::unique_ptr<RDKit::RWMol> d_mol;
  static constexpr unsigned int d_defaultWidth = 250;
  static constexpr unsigned int d_defaultHeight = 200;
};

class JSSubstructLibrary {
 public:
  JSSubstructLibrary(unsigned int num_bits);
  JSSubstructLibrary() : JSSubstructLibrary(d_defaultNumBits) {}
  int add_mol(const JSMol &m);
  int add_smiles(const std::string &smi);
  int add_trusted_smiles(const std::string &smi);
  JSMol *get_mol(unsigned int i);
  std::string get_matches(const JSMol &q, bool useChirality, int numThreads,
                          int maxResults) const;
  std::string get_matches(const JSMol &q, int maxResults) const {
    return get_matches(q, d_defaultUseChirality, d_defaultNumThreads,
                       maxResults);
  }
  std::string get_matches(const JSMol &q) const {
    return get_matches(q, d_defaultUseChirality, d_defaultNumThreads,
                       d_defaultMaxResults);
  }
  unsigned int count_matches(const JSMol &q, bool useChirality,
                             int numThreads) const;
  unsigned int count_matches(const JSMol &q) const {
    return count_matches(q, d_defaultUseChirality, d_defaultNumThreads);
  }

  std::unique_ptr<RDKit::SubstructLibrary> d_sslib;
  RDKit::CachedTrustedSmilesMolHolder *d_molHolder;
  RDKit::PatternHolder *d_fpHolder;
  unsigned int d_num_bits;
  static constexpr unsigned int d_defaultNumBits = 2048;
  static constexpr bool d_defaultUseChirality = true;
  static constexpr int d_defaultNumThreads = -1;
  static constexpr int d_defaultMaxResults = 1000;

 private:
  inline int add_mol_helper(const RDKit::ROMol &mol);
};

std::string get_inchikey_for_inchi(const std::string &input);
JSMol *get_mol(const std::string &input, const std::string &details_json);
JSMol *get_qmol(const std::string &input);
std::string version();
void prefer_coordgen(bool prefer);
