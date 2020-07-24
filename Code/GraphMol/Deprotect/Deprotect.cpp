//
//  Copyright (C) 2020 Brian P Kelley, Joann Prescott-Roy
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

namespace RDKit {
namespace Deprotect {
DeprotectData::DeprotectData(const std::string &deprotection_class,
                             const std::string &reaction_smarts,
                             const std::string &abbreviation,
                             const std::string &full_name)
    :

      deprotection_class(deprotection_class),
      reaction_smarts(reaction_smarts),
      abbreviation(abbreviation),
      full_name(full_name),
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
       "tert-butyldimethylsilyl"},
      {"alcohol", "COc1ccc(C([O;H0:1])(c2ccccc2)c2ccccc2)cc1>>[O;H1:1]", "MMT",
       "methoxytrityl"},
      {"alcohol", "O1C([O;H0:1])[C;H2][C;H2][C;H2][C;H2]1>>[O;H1:1]", "THP",
       "tetrahydropyranyl"},
      {"alcohol", "[C;H3][O;X2&R0][C;H2&R0][O;H0:1]>>[O;H1:1]", "MOM",
       "methoxymethyl_ether"},
      {"alcohol", "[C;R0][O;R0][C;R0][C;R0][O;R0][C;R0][O;X2:1]>>[O;H1:1]",
       "MEM", "β-Methoxyethoxymethyl_ether"},
      {"alcohol", "[O;!$(*C(=O)):1]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
       "Bn", "benzyl"},
      {"alcohol", "[O;H0&X2:1][Si]([C;H3])([C;H3])([C;H3])>>[O;H1:1]", "TMS",
       "trimethylsilyl"},
      {"alcohol", "[O;H0:1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[O;H1:1]",
       "Bz", "benzoyl"},
      {"alcohol", "[O;H0:1][C;R0](=O)[C;R0]([C;H3])([C;H3])[C;H3]>>[O;H1:1]",
       "Piv", "pivaloyl"},
      {"alcohol", "[O;H0:1][Si](C(C)C)(C(C)C)(C(C)C)>>[O;H1:1]", "TIPS",
       "triisopropylsilyl"},
      {"alcohol", "[O;R0:1][C;R0](=O)[C;H3]>>[O:1]", "Ac", "acetyl"},
      {"alcohol",
       "[c;H1]1[c;H1]c(O[C;H3])[c;H1][c;H1]c1[C;H2][O;D2&R0:1]>>[O;H1:1]",
       "PMB", "para-methoxybenzyl_ether"},
      {"alcohol", "c1ccc(C([NX3;H0,H1,H2:1])(c2ccccc2)c2ccccc2)cc1>>[N:1]",
       "Tr", "triphenylmethyl"},
      {"amine", "O=C1[N;H0:1]C(c2[c;H1][c;H1][c;H1][c;H1]c21)=O>>[N:1]", "Phth",
       "phthaloyl"},
      {"amine", "[#7:1]COCC[Si]([C;R0])([C;R0])[C;R0]>>[#7:1]", "SEM",
       "trimethylsilylethoxymethyl"},
      {"amine",
       "[C;H3;R0][O;R0]c1[c;H1][c;H1][c]([NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>"
       "[N:1]",
       "PMP", "para-methoxyphenyl"},
      {"amine",
       "[C;H3]c1[c;H1][c;H1]c(S(=O)(=O)[NX3;H0,H1;!$(NC=O):1])[c;H1][c;H1]1>>["
       "N:1]",
       "Ts", "tosyl"},
      {"amine",
       "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]",
       "Boc", "tert-butyloxycarbonyl"},
      {"amine", "[N;H0,H1:1]C(=O)C(F)(F)F>>[N:1]", "TFA", "trifluoroacetyl"},
      {"amine",
       "[NX3;H0,H1;!$(NC=O):1]C(=O)c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
       "Bz", "benzoyl"},
      {"amine", "[NX3;H0,H1;!$(NC=O):1]C(OCC1C(cccc2)c2c3c1cccc3)=O>>[N:1]",
       "Fmoc", "9-fluorenylmethyloxycarbonyl"},
      {"amine",
       "[NX3;H0,H1;!$(NC=O):1][C;H2]c1[c;H1][c;H1][c;H1][c;H1][c;H1]1>>[N:1]",
       "Bn", "benzyl"},
      {"amine", "[NX3;H0,H1;!$(NC=O):1][C;R0](=O)[C;H3]>>[N:1]", "Ac",
       "acetamide"},
      {"amine",
       "[NX3;H0,H1;!$(NC=O):1][C;R0](=O)[O;R0][C;R0]c1[c;H1][c;H1][c;H1][c;H1]["
       "c;H1]1>>[N:1]",
       "Cbz", "carbobenzyloxy"}};
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
          RDLog::BlockLogs blocker;
          MolOps::sanitizeMol(*dynamic_cast<RWMol *>(m.get()));
        } catch (MolSanitizeException &) {
          continue;
        }
        deprotections_used.push_back(deprotect.abbreviation);
        if (deprotections_used.size() >= MAX_DEPROTECTIONS) {
          BOOST_LOG(rdErrorLog)
              << "Too many deprotections, halting..." << std::endl;
        } else
          something_happened = true;
        break;
      }
    }
  }

  m->setProp("DEPROTECTIONS", deprotections_used);
  m->setProp<int>("DEPROTECTION_COUNT", deprotections_used.size());
  return std::unique_ptr<ROMol>(new ROMol(*m.get()));
};

}  // namespace Deprotect
} // namespace RDKit
