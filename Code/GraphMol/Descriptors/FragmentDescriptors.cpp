//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "FragmentDescriptors.h"

#include <memory>
#include <unordered_map>

namespace RDKit {
namespace Descriptors {
namespace {
struct FragmentDef {
  const char *name;
  const char *description;
  const char *smarts;
};

// generated from $RDBASE/Data/FragmentDescriptors.csv
const FragmentDef fragmentDefs[] = {
    {"fr_C_O", "Number of carbonyl O", "[CX3]=[OX1]"},
    {"fr_C_O_noCOO", "Number of carbonyl O, excluding COOH", "[C!$(C-[OH])]=O"},
    {"fr_Al_OH", "Number of aliphatic hydroxyl groups", "[C!$(C=O)]-[OH]"},
    {"fr_Ar_OH", "Number of aromatic hydroxyl groups", "c[OH1]"},
    {"fr_methoxy", "Number of methoxy groups -OCH3", "[OX2](-[#6])-[CH3]"},
    {"fr_oxime", "Number of oxime groups", "[CX3]=[NX2]-[OX2]"},
    {"fr_ester", "Number of esters", "[#6][CX3](=O)[OX2H0][#6]"},
    {"fr_Al_COO", "Number of aliphatic carboxylic acids", "C-C(=O)[O;H1,-]"},
    {"fr_Ar_COO", "Number of aromatic carboxylic acids", "c-C(=O)[O;H1,-]"},
    {"fr_COO", "Number of carboxylic acids", "[#6]C(=O)[O;H,-1]"},
    {"fr_COO2", "Number of carboxylic acids", "[CX3](=O)[OX1H0-,OX2H1]"},
    {"fr_ketone", "Number of ketones", "[#6][CX3](=O)[#6]"},
    {"fr_ether", "Number of ether oxygens (including phenoxy)",
     "[OD2]([#6])[#6]"},
    {"fr_phenol", "Number of phenols", "[OX2H]-c1ccccc1"},
    {"fr_aldehyde", "Number of aldehydes", "[CX3H1](=O)[#6]"},
    {"fr_quatN", "Number of quaternary nitrogens", "[$([NX4+]),$([NX4]=*)]"},
    {"fr_NH2", "Number of Primary amines", "[NH2,nH2]"},
    {"fr_NH1", "Number of Secondary amines", "[NH1,nH1]"},
    {"fr_NH0", "Number of Tertiary amines", "[NH0,nH0]"},
    {"fr_Ar_N", "Number of aromatic nitrogens", "n"},
    {"fr_Ar_NH", "Number of aromatic amines", "[nH]"},
    {"fr_aniline", "Number of anilines", "c-[NX3;!$(N=*)]"},
    {"fr_Imine", "Number of Imines", "[Nv3](=C)-[#6]"},
    {"fr_nitrile", "Number of nitriles", "[NX1]#[CX2]"},
    {"fr_hdrzine", "Number of hydrazine groups", "[NX3]-[NX3]"},
    {"fr_hdrzone", "Number of hydrazone groups", "C=N-[NX3]"},
    {"fr_nitroso", "Number of nitroso groups, excluding NO2", "[N!$(N-O)]=O"},
    {"fr_N_O", "Number of hydroxylamine groups", "[N!$(N=O)](-O)-C"},
    {"fr_nitro", "Number of nitro groups",
     "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"},
    {"fr_azo", "Number of azo groups", "[#6]-N=N-[#6]"},
    {"fr_diazo", "Number of diazo groups", "[N+]#N"},
    {"fr_azide", "Number of azide groups",
     "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]"},
    {"fr_amide", "Number of amides", "C(=O)-N"},
    {"fr_priamide", "Number of primary amides", "C(=O)-[NH2]"},
    {"fr_amidine", "Number of amidine groups", "C(=N)(-N)-[!#7]"},
    {"fr_guanido", "Number of guanidine groups", "C(=N)(N)N"},
    {"fr_Nhpyrrole", "Number of H-pyrrole nitrogens", "[nH]"},
    {"fr_imide", "Number of imide groups", "N(-C(=O))-C=O"},
    {"fr_isocyan", "Number of isocyanates", "N=C=O"},
    {"fr_isothiocyan", "Number of isothiocyanates", "N=C=S"},
    {"fr_thiocyan", "Number of thiocyanates", "S-C#N"},
    {"fr_halogen", "Number of halogens", "[#9,#17,#35,#53]"},
    {"fr_alkyl_halide", "Number of alkyl halides", "[CX4]-[Cl,Br,I,F]"},
    {"fr_sulfide", "Number of thioether", "[SX2](-[#6])-C"},
    {"fr_SH", "Number of thiol groups", "[SH]"},
    {"fr_C_S", "Number of thiocarbonyl", "C=[SX1]"},
    {"fr_sulfone", "Number of sulfone groups",
     "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]"},
    {"fr_sulfonamd", "Number of sulfonamides",
     "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]"},
    {"fr_prisulfonamd", "Number of primary sulfonamides",
     "[NH2]-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]"},
    {"fr_barbitur", "Number of barbiturate groups", "C1C(=O)NC(=O)NC1=O"},
    {"fr_urea", "Number of urea groups", "C(=O)(-N)-N"},
    {"fr_term_acetylene", "Number of terminal acetylenes", "C#[CH]"},
    {"fr_imidazole", "Number of imidazole rings", "n1cncc1"},
    {"fr_furan", "Number of furan rings", "o1cccc1"},
    {"fr_thiophene", "Number of thiophene rings", "s1cccc1"},
    {"fr_thiazole", "Number of thiazole rings", "c1scnc1"},
    {"fr_oxazole", "Number of oxazole rings", "c1ocnc1"},
    {"fr_pyridine", "Number of pyridine rings", "n1ccccc1"},
    {"fr_piperdine", "Number of piperdine rings", "N1CCCCC1"},
    {"fr_piperzine", "Number of piperzine rings", "N1CCNCC1"},
    {"fr_morpholine", "Number of morpholine rings", "O1CCNCC1"},
    {"fr_lactam", "Number of beta lactams", "N1C(=O)CC1"},
    {"fr_lactone", "Number of cyclic esters (lactones)",
     "[C&R1](=O)[O&R1][C&R1]"},
    {"fr_tetrazole", "Number of tetrazole rings", "c1nnnn1"},
    {"fr_epoxide", "Number of epoxide rings", "O1CC1"},
    {"fr_unbrch_alkane",
     "Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)",
     "[CR0;D2,D1][CR0;D2][CR0;D2][CR0;D2,D1]"},
    {"fr_bicyclic", "Bicyclic", "[R2][R2]"},
    {"fr_benzene", "Number of benzene rings", "c1ccccc1"},
    {"fr_phos_acid", "Number of phosphoric acid groups",
     "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"},
    {"fr_phos_ester", "Number of phosphoric ester groups",
     "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"},
    {"fr_nitro_arom", "Number of nitro benzene ring substituents",
     "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1)]"},
    {"fr_nitro_arom_nonortho",
     "Number of non-ortho nitro benzene ring substituents",
     "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]"},
    {"fr_dihydropyridine", "Number of dihydropyridines",
     "[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]"},
    {"fr_phenol_noOrthoHbond",
     "Number of phenolic OH excluding ortho intramolecular Hbond substituents",
     "[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]"},
    {"fr_Al_OH_noTert", "Number of aliphatic hydroxyl groups excluding tert-OH",
     "[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]"},
    {"fr_benzodiazepine",
     "Number of benzodiazepines with no additional fused rings",
     "[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2"},
    {"fr_para_hydroxylation", "Number of para-hydroxylation sites",
     "[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]"},
    {"fr_allylic_oxid",
     "Number of allylic oxidation sites excluding steroid dienone",
     "[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]"},
    {"fr_aryl_methyl", "Number of aryl methyl sites for hydroxylation",
     "[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]"},
    {"fr_Ndealkylation1", "Number of XCCNR groups",
     "[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]"},
    {"fr_Ndealkylation2",
     "Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)",
     "[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]"},
    {"fr_alkyl_carbamate", "Number of alkyl carbamates (subject to hydrolysis)",
     "C[NH1]C(=O)OC"},
    {"fr_ketone_Topliss",
     "Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha",
     "[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]"},
    {"fr_ArN", "Number of N functional groups attached to aromatics",
     "[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]"},
    {"fr_HOCCN", "Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic",
     "[$([OX2H1][CX4][CX4H2][NX3&R1]),$([OH1][CX4][CX4H2][NX3][CX4](C)(C)C)]"},
};

struct FragmentPatterns {
  std::vector<std::string> names;
  std::vector<std::string> descriptions;
  std::vector<std::string> smarts;
  std::vector<std::unique_ptr<ROMol>> patterns;
  std::unordered_map<std::string, size_t> nameIdx;
  FragmentPatterns() {
    for (const auto &def : fragmentDefs) {
      std::unique_ptr<ROMol> patt(SmartsToMol(def.smarts));
      CHECK_INVARIANT(
          patt && patt->getNumAtoms(),
          std::string("could not parse fragment SMARTS ") + def.smarts);
      nameIdx[def.name] = names.size();
      names.emplace_back(def.name);
      descriptions.emplace_back(def.description);
      smarts.emplace_back(def.smarts);
      patterns.push_back(std::move(patt));
    }
  }
};

const FragmentPatterns &getPatterns() {
  static const FragmentPatterns patterns;
  return patterns;
}

unsigned int countMatches(const ROMol &mol, const ROMol &patt,
                          bool countUnique) {
  SubstructMatchParameters ps;
  ps.uniquify = countUnique;
  return SubstructMatch(mol, patt, ps).size();
}
}  // namespace

const std::vector<std::string> &getFragmentDescriptorNames() {
  return getPatterns().names;
}

const std::vector<std::string> &getFragmentDescriptorSmarts() {
  return getPatterns().smarts;
}

const std::vector<std::string> &getFragmentDescriptorDescriptions() {
  return getPatterns().descriptions;
}

std::vector<unsigned int> calcFragmentDescriptors(const ROMol &mol,
                                                  bool countUnique) {
  const auto &patts = getPatterns();
  std::vector<unsigned int> res;
  res.reserve(patts.patterns.size());
  for (const auto &patt : patts.patterns) {
    res.push_back(countMatches(mol, *patt, countUnique));
  }
  return res;
}

unsigned int calcFragmentDescriptor(const ROMol &mol, const std::string &name,
                                    bool countUnique) {
  const auto &patts = getPatterns();
  auto it = patts.nameIdx.find(name);
  if (it == patts.nameIdx.end()) {
    throw KeyErrorException(name);
  }
  return countMatches(mol, *patts.patterns[it->second], countUnique);
}

}  // namespace Descriptors
}  // namespace RDKit
