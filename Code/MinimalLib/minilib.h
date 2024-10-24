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
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#ifdef RDK_BUILD_MINIMAL_LIB_MMPA
#include <GraphMol/MMPA/MMPA.h>
#endif
#ifdef RDK_BUILD_MINIMAL_LIB_RGROUPDECOMP
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#endif

class JSMolList;

class JSMolBase {
 public:
  JSMolBase(const JSMolBase &) = delete;
  JSMolBase &operator=(const JSMolBase &) = delete;
  virtual ~JSMolBase() {};
  virtual const RDKit::RWMol &get() const = 0;
  virtual RDKit::RWMol &get() = 0;
  std::string get_smiles() const;
  std::string get_smiles(const std::string &details) const;
  std::string get_cxsmiles() const;
  std::string get_cxsmiles(const std::string &details) const;
  std::string get_smarts() const;
  std::string get_smarts(const std::string &details) const;
  std::string get_cxsmarts() const;
  std::string get_cxsmarts(const std::string &details) const;
  std::string get_molblock(const std::string &details) const;
  std::string get_molblock() const { return get_molblock("{}"); }
  std::string get_v3Kmolblock(const std::string &details) const;
  std::string get_v3Kmolblock() const { return get_v3Kmolblock("{}"); }
  std::string get_pickle() const;
#ifdef RDK_BUILD_INCHI_SUPPORT
  std::string get_inchi(const std::string &options) const;
  std::string get_inchi() const { return get_inchi(""); }
#endif
  std::string get_json() const;
  std::string get_svg(int width, int height) const;
  std::string get_svg() const {
    return get_svg(d_defaultWidth, d_defaultHeight);
  }
  std::string get_svg_with_highlights(const std::string &details) const;
  std::string get_substruct_match(const JSMolBase &q) const;
  std::string get_substruct_matches(const JSMolBase &q) const;
  std::string get_descriptors() const;
  std::string get_morgan_fp(const std::string &details) const;
  std::string get_morgan_fp() const { return get_morgan_fp("{}"); }
  std::string get_morgan_fp_as_binary_text(const std::string &details) const;
  std::string get_morgan_fp_as_binary_text() const {
    return get_morgan_fp_as_binary_text("{}");
  }
  std::string get_pattern_fp(const std::string &details) const;
  std::string get_pattern_fp() const { return get_pattern_fp("{}"); }
  std::string get_pattern_fp_as_binary_text(const std::string &details) const;
  std::string get_pattern_fp_as_binary_text() const {
    return get_pattern_fp_as_binary_text("{}");
  }
  std::string get_topological_torsion_fp(const std::string &details) const;
  std::string get_topological_torsion_fp() const {
    return get_topological_torsion_fp("{}");
  };
  std::string get_topological_torsion_fp_as_binary_text(
      const std::string &details) const;
  std::string get_topological_torsion_fp_as_binary_text() const {
    return get_topological_torsion_fp_as_binary_text("{}");
  }
  std::string get_rdkit_fp(const std::string &details) const;
  std::string get_rdkit_fp() const { return get_rdkit_fp("{}"); }
  std::string get_rdkit_fp_as_binary_text(const std::string &details) const;
  std::string get_rdkit_fp_as_binary_text() const {
    return get_rdkit_fp_as_binary_text("{}");
  }
  std::string get_atom_pair_fp(const std::string &details) const;
  std::string get_atom_pair_fp() const { return get_atom_pair_fp("{}"); }
  std::string get_atom_pair_fp_as_binary_text(const std::string &details) const;
  std::string get_atom_pair_fp_as_binary_text() const {
    return get_atom_pair_fp_as_binary_text("{}");
  }
  std::string get_maccs_fp() const;
  std::string get_maccs_fp_as_binary_text() const;
#ifdef RDK_BUILD_AVALON_SUPPORT
  std::string get_avalon_fp(const std::string &details) const;
  std::string get_avalon_fp() const { return get_avalon_fp("{}"); }
  std::string get_avalon_fp_as_binary_text(const std::string &details) const;
  std::string get_avalon_fp_as_binary_text() const {
    return get_avalon_fp_as_binary_text("{}");
  }
#endif
  std::string condense_abbreviations(double maxCoverage, bool useLinkers);
  std::string condense_abbreviations() {
    return condense_abbreviations(0.4, false);
  }
  std::string condense_abbreviations_from_defs(const std::string &definitions,
                                               double maxCoverage,
                                               bool areLinkers);
  std::string generate_aligned_coords(const JSMolBase &templateMol,
                                      const std::string &details);
  std::string generate_aligned_coords(const JSMolBase &templateMol) {
    return generate_aligned_coords(templateMol, "{}");
  }
  [[deprecated(
      "please check the get_mol/get_qmol return value for non-nullness "
      "instead")]] bool
  is_valid() const;
  int has_coords() const;

  std::string get_stereo_tags();
  std::string get_aromatic_form() const;
  void convert_to_aromatic_form();
  std::string get_kekule_form() const;
  void convert_to_kekule_form();
  bool set_new_coords(bool useCoordGen);
  bool set_new_coords() { return set_new_coords(false); }
  std::string get_new_coords(bool useCoordGen) const;
  std::string get_new_coords() const { return get_new_coords(false); }
  bool has_prop(const std::string &key) const;
  std::vector<std::string> get_prop_list(bool includePrivate,
                                         bool includeComputed) const;
  std::vector<std::string> get_prop_list(bool includePrivate) const {
    return get_prop_list(includePrivate, true);
  }
  std::vector<std::string> get_prop_list() const {
    return get_prop_list(true, true);
  }
  bool set_prop(const std::string &key, const std::string &val, bool computed);
  bool set_prop(const std::string &key, const std::string &val) {
    return set_prop(key, val, false);
  }
  std::string get_prop(const std::string &key) const;
  bool clear_prop(const std::string &key);
  std::string remove_hs() const;
  bool remove_hs_in_place();
  std::string add_hs() const;
  bool add_hs_in_place();
  double normalize_depiction(int canonicalize, double scaleFactor);
  double normalize_depiction(int canonicalize) {
    return normalize_depiction(canonicalize, -1.);
  }
  double normalize_depiction() { return normalize_depiction(1, -1.); }
  void straighten_depiction(bool minimizeRotation);
  void straighten_depiction() { straighten_depiction(false); }
  std::pair<JSMolList *, std::string> get_frags(
      const std::string &details_json) const;
  std::pair<JSMolList *, std::string> get_frags() const {
    return get_frags("{}");
  }
  unsigned int get_num_atoms(bool heavyOnly) const;
  unsigned int get_num_atoms() const { return get_num_atoms(false); };
  unsigned int get_num_bonds() const;
#ifdef RDK_BUILD_MINIMAL_LIB_MMPA
  std::pair<JSMolList *, JSMolList *> get_mmpa_frags(
      unsigned int minCuts, unsigned int maxCuts,
      unsigned int maxCutBonds) const;
#endif

  static constexpr int d_defaultWidth = 250;
  static constexpr int d_defaultHeight = 200;

 protected:
  JSMolBase() = default;
};

class JSMol : public JSMolBase {
 public:
  JSMol() : d_mol(new RDKit::RWMol()) {}
  JSMol(RDKit::RWMol *mol) : d_mol(mol) { checkNotNull(); }
  JSMol(const JSMol &other) { d_mol.reset(new RDKit::RWMol(other.get())); }
  JSMol &operator=(const JSMol &other) {
    d_mol.reset(new RDKit::RWMol(other.get()));
    return *this;
  }
  const RDKit::RWMol &get() const {
    checkNotNull();
    return *d_mol.get();
  }
  RDKit::RWMol &get() {
    checkNotNull();
    return *d_mol.get();
  }

 private:
  void checkNotNull() const { CHECK_INVARIANT(d_mol, "d_mol cannot be null"); }
  std::unique_ptr<RDKit::RWMol> d_mol;
};

class JSMolShared : public JSMolBase {
 public:
  JSMolShared() = delete;
  JSMolShared(const RDKit::ROMOL_SPTR &mol) : d_mol(mol) { checkNotNull(); }
  JSMolShared(const JSMolShared &other) { d_mol = other.d_mol; }
  JSMolShared &operator=(const JSMolShared &other) {
    d_mol = other.d_mol;
    return *this;
  }
  const RDKit::RWMol &get() const {
    checkNotNull();
    return static_cast<RDKit::RWMol &>(*d_mol.get());
  }
  RDKit::RWMol &get() {
    checkNotNull();
    return static_cast<RDKit::RWMol &>(*d_mol.get());
  }
  const RDKit::ROMOL_SPTR &get_sptr() const { return d_mol; }
  RDKit::ROMOL_SPTR &get_sptr() { return d_mol; }

 private:
  void checkNotNull() const { CHECK_INVARIANT(d_mol, "d_mol cannot be null"); }
  RDKit::ROMOL_SPTR d_mol;
};

class JSMolList {
 public:
  JSMolList(const std::vector<RDKit::ROMOL_SPTR> &mols)
      : d_mols(mols), d_idx(0) {};
  JSMolList() : d_idx(0) {};
  JSMolBase *next();
  size_t append(const JSMolBase &mol);
  size_t insert(size_t idx, const JSMolBase &mol);
  JSMolBase *at(size_t idx) const;
  JSMolBase *pop(size_t idx);
  void reset() { d_idx = 0; }
  bool at_end() const { return d_idx == d_mols.size(); }
  size_t size() const { return d_mols.size(); }
  const std::vector<RDKit::ROMOL_SPTR> &mols() const { return d_mols; }

 private:
  std::vector<RDKit::ROMOL_SPTR> d_mols;
  size_t d_idx;
};

namespace RDKit {
namespace MinimalLib {
struct LogHandle;
}
}  // namespace RDKit

class JSLog {
 public:
  JSLog(RDKit::MinimalLib::LogHandle *logHandle);
  ~JSLog();
  std::string get_buffer() const;
  void clear_buffer() const;

 private:
  RDKit::MinimalLib::LogHandle *d_logHandle;
};

#ifdef RDK_BUILD_MINIMAL_LIB_RXN
class JSReaction {
 public:
  JSReaction() : d_rxn(new RDKit::ChemicalReaction()) {}
  JSReaction(RDKit::ChemicalReaction *rxn) : d_rxn(rxn) { assert(d_rxn); }
  [[deprecated(
      "please check the get_rxn return value for non-nullness "
      "instead")]] bool
  is_valid() const;

  std::vector<JSMolList *> run_reactants(const JSMolList &reactants,
                                         unsigned int maxProducts) const;
  static constexpr int maxProducts = 1000;
  std::string get_svg(int width, int height) const;
  std::string get_svg() const {
    return get_svg(d_defaultWidth, d_defaultHeight);
  }
  std::string get_svg_with_highlights(const std::string &details) const;

  std::unique_ptr<RDKit::ChemicalReaction> d_rxn;
  static constexpr int d_defaultWidth = 800;
  static constexpr int d_defaultHeight = 200;
};
#endif

#ifdef RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY
class JSSubstructLibrary {
 public:
  JSSubstructLibrary(unsigned int num_bits);
  JSSubstructLibrary() : JSSubstructLibrary(d_defaultNumBits) {}
  int add_mol(const JSMolBase &m);
  int add_smiles(const std::string &smi);
  int add_trusted_smiles(const std::string &smi);
  int add_trusted_smiles_and_pattern_fp(const std::string &smi,
                                        const std::string &patternFp);
  std::string get_trusted_smiles(unsigned int i) const;
  std::string get_pattern_fp(unsigned int i) const;
  JSMolBase *get_mol(unsigned int i);
  std::string get_matches(const JSMolBase &q, bool useChirality, int numThreads,
                          int maxResults) const;
  std::string get_matches(const JSMolBase &q, int maxResults) const {
    return get_matches(q, d_defaultUseChirality, d_defaultNumThreads,
                       maxResults);
  }
  std::string get_matches(const JSMolBase &q) const {
    return get_matches(q, d_defaultUseChirality, d_defaultNumThreads,
                       d_defaultMaxResults);
  }
  unsigned int count_matches(const JSMolBase &q, bool useChirality,
                             int numThreads) const;
  unsigned int count_matches(const JSMolBase &q, bool useChirality) const {
    return count_matches(q, useChirality, d_defaultNumThreads);
  }
  unsigned int count_matches(const JSMolBase &q) const {
    return count_matches(q, d_defaultUseChirality, d_defaultNumThreads);
  }
  unsigned int size() const { return d_sslib->size(); }

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
#endif

#ifdef RDK_BUILD_INCHI_SUPPORT
std::string get_inchikey_for_inchi(const std::string &input);
#endif
JSMolBase *get_mol(const std::string &input, const std::string &details_json);
JSMolBase *get_mol_from_pickle(const std::string &pkl);
JSMolBase *get_mol_copy(const JSMolBase &other);
JSMolBase *get_qmol(const std::string &input);
#ifdef RDK_BUILD_MINIMAL_LIB_RXN
JSReaction *get_rxn(const std::string &input, const std::string &details_json);
#endif
std::string version();
void prefer_coordgen(bool prefer);
bool use_legacy_stereo_perception(bool value);
bool allow_non_tetrahedral_chirality(bool value);
void enable_logging();
void disable_logging();
JSLog *set_log_tee(const std::string &log_name);
JSLog *set_log_capture(const std::string &log_name);
#ifdef RDK_BUILD_MINIMAL_LIB_MCS
std::string get_mcs_as_json(const JSMolList &mols,
                            const std::string &details_json);
std::string get_mcs_as_smarts(const JSMolList &mols,
                              const std::string &details_json);
JSMolBase *get_mcs_as_mol(const JSMolList &mols,
                          const std::string &details_json);
#endif

#ifdef RDK_BUILD_MINIMAL_LIB_RGROUPDECOMP
class JSRGroupDecomposition {
 public:
  JSRGroupDecomposition(const JSMolBase &core, const std::string &details_json);
  JSRGroupDecomposition(const JSMolBase &core)
      : JSRGroupDecomposition(core, "") {};
  JSRGroupDecomposition(const JSMolList &cores,
                        const std::string &details_json);
  JSRGroupDecomposition(const JSMolList &cores)
      : JSRGroupDecomposition(cores, "") {};
  int add(const JSMolBase &mol);
  bool process();
  std::map<std::string, std::unique_ptr<JSMolList>> getRGroupsAsColumns() const;
  std::vector<std::map<std::string, std::unique_ptr<JSMolBase>>>
  getRGroupsAsRows() const;

 private:
  std::unique_ptr<RDKit::RGroupDecomposition> d_decomp;
  std::vector<unsigned int> d_unmatched;
};
#endif
