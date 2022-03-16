//
//  Copyright (C) 2020-2021 Brian P Kelley, Joann Prescott-Roy and other RDKit
//  contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Deprotect.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <boost/smart_ptr.hpp>
#include <utility>

namespace RDKit {
namespace Deprotect {
DeprotectData::DeprotectData(std::string deprotection_class,
                             const std::string &reaction_smarts,
                             std::string abbreviation, std::string full_name,
                             std::string example)
    : deprotection_class(std::move(deprotection_class)),
      reaction_smarts(reaction_smarts),
      abbreviation(std::move(abbreviation)),
      full_name(std::move(full_name)),
      example(std::move(example)),
      rxn(RxnSmartsToChemicalReaction(reaction_smarts)) {
  if (rxn.get()) {
    if (rxn->getNumProductTemplates() != 1) {
      BOOST_LOG(rdErrorLog)
          << "Deprotection reactions must have exactly one product"
          << std::endl;
    }
    rxn->initReactantMatchers();
  }
}

const std::vector<DeprotectData> &getDeprotections() {
  static const std::vector<DeprotectData> deprotections{
      {"alcohol", "CC(C)([Si](C)(C)[O;H0:1])C>>[O;H1:1]", "TBDMS",
       "tert-butyldimethylsilyl", "CC(C)(C)[Si](C)(C)OC>>OC"},
      {"alcohol", "COc1ccc(C([O;H0:1])(c2ccccc2)c2ccccc2)cc1>>[O;H1:1]", "MMT",
       "methoxytrityl", "COc1ccc(C(ON)(c2ccccc2)c2ccccc2)cc1>>NO"},
      {"alcohol", "O1C([O;H0:1])[C;H2][C;H2][C;H2][C;H2]1>>[O;H1:1]", "THP",
       "tetrahydropyranyl", "COC1CCCCO1>>CO"},
      {"alcohol", "[C;H3][O;X2&R0][C;H2&R0][O;H0:1]>>[O;H1:1]", "MOM",
       "methoxymethyl_ether", "COCOC>>CO"},
      {"alcohol", "[C;R0][O;R0][C;R0][C;R0][O;R0][C;R0][O;X2:1]>>[O;H1:1]",
       "MEM", "beta-Methoxyethoxymethyl_ether", "C(C)(C)OCOCCOC>>C(C)(C)O"},
      {"alcohol", "[O;!$(*C(=O)):1]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
       "Bn", "benzyl", "NOc1ccccc1>>ON"},
      {"alcohol", "[O;H0&X2:1][Si]([C;H3])([C;H3])([C;H3])>>[O;H1:1]", "TMS",
       "trimethylsilyl", "NO[Si](C)(C)C>>NO"},
      {"alcohol", "[O;H0:1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
       "Bz", "benzoyl", "COC(=O)c1ccccc1>>OC"},
      {"alcohol", "[O;H0:1][C;R0](=O)[C;R0]([C;H3])([C;H3])[C;H3]>>[O;H1:1]",
       "Piv", "pivaloyl", "CC(C)(C)C(=O)ON>>NO"},
      {"alcohol", "[O;H0:1][Si](C(C)C)(C(C)C)(C(C)C)>>[O;H1:1]", "TIPS",
       "triisopropylsilyl", "NO[Si](C(C)C)(C(C)C)(C(C)C)>>ON"},
      {"alcohol", "[O;R0:1][C;R0](=O)[C;H3]>>[O:1]", "Ac", "acetyl",
       "NOC(=O)C>>NO"},
      {"alcohol",
       "[c;H1]1[c;H1]c(O[C;H3])[c;H1][c;H1]c1[C;H2][O;D2&R0:1]>>[O;H1:1]",
       "PMB", "para-methoxybenzyl_ether", "c1cc(OC)ccc1CON>>NO"},
      {"alcohol", "c1ccc(C([NX3;H0,H1,H2:1])(c2ccccc2)c2ccccc2)cc1>>[N:1]",
       "Tr", "triphenylmethyl", "CN(C)C(c1ccccc1)(c1ccccc1)c1ccccc1>>CNC"},
      {"amine", "O=C1[N;H0:1]C(c2[c;H1][c;H1][c;H1][c;H1]c21)=O>>[N:1]", "Phth",
       "phthaloyl", "O=C1N(C(=O)c2ccccc12)C1=CC=CC=C1>>Nc1ccccc1"},
      {"amine", "[#7:1]COCC[Si]([C;R0])([C;R0])[C;R0]>>[#7:1]", "SEM",
       "trimethylsilylethoxymethyl", "C[Si](C)(C)CCOCNO>>NO"},
      {"amine",
       "[C;H3;R0][O;R0]c1[c;H1][c;H1][c]([NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>"
       "[N:1]",
       "PMP", "para-methoxyphenyl", "COc1ccc(NO)cc1>>NO"},
      {"amine",
       "[C;H3]c1[c;H1][c;H1]c(S(=O)(=O)[NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>["
       "N:1]",
       "Ts", "tosyl", "Cc1ccc(cc1)S(NO)(=O)=O>>NO"},
      {"amine",
       "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]",
       "Boc", "tert-butyloxycarbonyl",
       "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C>>NCc1ccccc1N"},
      {"amine", "[N;H0,H1:1]C(=O)C(F)(F)F>>[N:1]", "TFA", "trifluoroacetyl",
       "ONC(=O)C(F)(F)F>>NO"},
      {"amine", "[NX3;H0,H1:1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
       "Bz", "benzoyl", "ONC(=O)c1ccccc1>>NO"},
      {"amine",
       "[NX3;H0,H1:1][#6](=O)-[#8]-[#6]-[#6]-1-c2ccccc2-c2ccccc-12>>[N:1]",
       "Fmoc", "9-fluorenylmethyloxycarbonyl",
       "CN(C)C(=O)OCC1C2=CC=CC=C2C2=C1C=CC=C2>>CNC"},
      {"amine",
       "[NX3;H0,H1;!$(NC=O):1][C;H2]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
       "Bn", "benzylamine", "CN(C)Cc1ccccc1>>CNC"},
      {"amine", "[NX3;H0,H1:1][C;R0](=O)[C;H3]>>[N:1]", "Ac", "acetamide",
       "N(C)(C)C(=O)C>>CNC"},
      {"amine",
       "[NX3;H0,H1:1][C;R0](=O)[O;R0][C;R0]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>["
       "N:1]",
       "Cbz", "carbobenzyloxy", "N(C)(C)C(=O)OCc1ccccc1>>CNC"},
      {"carbonyl",
       "[#6,#1:1][C:2]([#6,#1:3])([#8:4][#6])([#8][#6])>>[*:1][C:2]([*:3])[#8:"
       "4]",
       "Acetyl/Ketal", "Acetal/Ketal",
       "O=C(OC)C1CCC(OC2)(OC2)CC1>>COC(=O)C1CCC(O)CC1"},
  };
  return deprotections;
}

const unsigned int MAX_DEPROTECTIONS = 100;
std::unique_ptr<ROMol> deprotect(
    const ROMol &mol, const std::vector<DeprotectData> &deprotections) {
  std::vector<std::string> deprotections_used;

  bool something_happened = true;
  auto m = boost::make_shared<ROMol>(mol);
  while (something_happened) {
    something_happened = false;
    for (auto &deprotect : deprotections) {
      if (!deprotect.isValid()) {
        // error and contine;
        continue;
      }
      for (auto &prods : deprotect.rxn->runReactant(m, 0)) {
        m = prods[0];
        try {
          RDLog::LogStateSetter blocker;
          MolOps::sanitizeMol(*dynamic_cast<RWMol *>(m.get()));
        } catch (MolSanitizeException &) {
          continue;
        }
        deprotections_used.push_back(deprotect.abbreviation);
        if (deprotections_used.size() >= MAX_DEPROTECTIONS) {
          BOOST_LOG(rdErrorLog)
              << "Too many deprotections, halting..." << std::endl;
        } else {
          something_happened = true;
        }
        break;
      }
    }
  }

  m->setProp("DEPROTECTIONS", deprotections_used);
  m->setProp<int>("DEPROTECTION_COUNT", deprotections_used.size());
  return std::unique_ptr<ROMol>(new ROMol(*m.get()));
}

bool deprotectInPlace(RWMol &mol,
                      const std::vector<DeprotectData> &deprotections) {
  std::vector<std::string> deprotections_used;

  bool modified = false;
  bool something_happened = true;
  while (something_happened) {
    something_happened = false;
    for (auto &deprotect : deprotections) {
      if (!deprotect.isValid()) {
        // error and contine;
        continue;
      }
      bool changes = deprotect.rxn->runReactant(mol);
      if (changes) {
        try {
          RDLog::LogStateSetter blocker;
          MolOps::sanitizeMol(mol);
        } catch (MolSanitizeException &) {
          continue;
        }
        modified = true;
        deprotections_used.push_back(deprotect.abbreviation);
        if (deprotections_used.size() >= MAX_DEPROTECTIONS) {
          BOOST_LOG(rdErrorLog)
              << "Too many deprotections, halting..." << std::endl;
        } else {
          something_happened = true;
        }
      }
    }
  }

  mol.setProp("DEPROTECTIONS", deprotections_used);
  mol.setProp<int>("DEPROTECTION_COUNT", deprotections_used.size());
  return modified;
};

}  // namespace Deprotect
}  // namespace RDKit
