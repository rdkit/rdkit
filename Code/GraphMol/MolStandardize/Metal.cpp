//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Metal.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>

using namespace std;
using namespace RDKit;
namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

MetalDisconnector::MetalDisconnector()
    : metal_nof(
          SmartsToMol("[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,"
                      "Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,"
                      "W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]")),
      metal_non(SmartsToMol(
          "[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,"
          "W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]")){
			BOOST_LOG(rdInfoLog) << "Initializing MetalDisconnector\n";
			};

MetalDisconnector::MetalDisconnector(const MetalDisconnector &other) {
  metal_nof = other.metal_nof;
  metal_non = other.metal_non;
};

MetalDisconnector::~MetalDisconnector(){};

ROMol *MetalDisconnector::getMetalNof() {return metal_nof.get();}

ROMol *MetalDisconnector::getMetalNon() {return metal_non.get();}

void MetalDisconnector::setMetalNof(const ROMol &mol) {
	this->metal_nof = static_cast<ROMOL_SPTR>(new ROMol(mol));
}

void MetalDisconnector::setMetalNon(const ROMol &mol) {
	this->metal_non = static_cast<ROMOL_SPTR>(new ROMol(mol));
}

ROMol *MetalDisconnector::disconnect(const ROMol &mol) {
  auto *res = new RWMol(mol);
  MetalDisconnector::disconnect(*res);
  return static_cast<ROMol *>(res);
}

void MetalDisconnector::disconnect(RWMol &mol) {
	BOOST_LOG(rdInfoLog) << "Running MetalDisconnector\n";
  std::list<ROMOL_SPTR> metalList = {metal_nof, metal_non};
  for (auto &query : metalList) {
    std::vector<MatchVectType> matches;
    unsigned int matched;

    matched = SubstructMatch(mol, *query, matches);

    for (size_t i = 0; i < matched; ++i) {
      // disconnecting metal-R bond
      int metal_idx = matches[i][0].second;
      int non_idx = matches[i][1].second;
      Bond *b = mol.getBondBetweenAtoms(metal_idx, non_idx);
      double order = b->getBondTypeAsDouble();
      mol.removeBond(metal_idx, non_idx);

      // adjusting neighbouring charges
      Atom *a1 = mol.getAtomWithIdx(metal_idx);
      a1->setFormalCharge(a1->getFormalCharge() + int(order));
      Atom *a2 = mol.getAtomWithIdx(non_idx);
      a2->setFormalCharge(a2->getFormalCharge() - int(order));
			BOOST_LOG(rdInfoLog) << "Removed covalent bond between "
				<< a1->getSymbol() << " and " << a2->getSymbol() << "\n";				
    }
    //	std::cout << "After removing bond and charge adjustment: " <<
    //MolToSmiles(mol) << std::endl;
  }
  MolOps::sanitizeMol(mol);
}
}  // namespace MolStandardize
}  // namespace RDKit
