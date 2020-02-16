//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_HAS_EIGEN3
#include "BCUT.h"
#include <Eigen/Dense>
#include <GraphMol/RDKitBase.h>
#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <RDGeneral/types.h>

namespace RDKit {
namespace Descriptors {
  // diagonal elements are a property (atomic num, charge, etc)
  // off diagonal are .1,.2,.3.15 for single,double,triple or aromatic
  // all other elements are .001
  std::pair<double,double> BCUT2D(const ROMol &m, const std::string &atom_double_prop) {
    unsigned int num_atoms = m.getNumAtoms();
    if (num_atoms == 0) {
        return std::pair<double,double>(0,0);
    }
    
    Eigen::MatrixXd burden(num_atoms, num_atoms);

    // All non bonded entires set to 0.001
    for(unsigned int i=0;i<num_atoms;++i) {
      burden(i,i) = m.getAtomWithIdx(i)->getProp<double>(atom_double_prop);
      for(unsigned int j=0;j<num_atoms;++j) {
	if (i!=j)
	  burden(i,j) = 0.001;
      }
    }

 for(auto &bond : m.bonds()) {
      unsigned int i = bond->getBeginAtomIdx();
      unsigned int j = bond->getEndAtomIdx();
      double score = 0.0;
      switch(bond->getBondType()) {
      case Bond::AROMATIC:
	// score = 0.15; orig burden
	score = 0.8164965809277261; // 1/sqrt(1.5)
	break;
      case Bond::SINGLE:
	// score = 0.1;
	score = 1.0; // 1/sqrt(1.0)
	break;
      case Bond::DOUBLE:
	// score = 0.2;
	score = 0.7071067811865475; // 1/sqrt(2.0)
	break;
      case Bond::TRIPLE:
	// score = 0.3;	
	score = 0.5773502691896258; // 1/sqrt(3);
	break;
      default:
	CHECK_INVARIANT(0, "Bond order must be Single, Double, Triple or Aromatic");
      }
      burden(i,j) = burden(j,i) = score;
    }    
    
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    std::cerr << burden.format(OctaveFmt) << std::endl;
    Eigen::VectorXcd eivals = burden.eigenvalues();
    double highest, lowest;
    highest = lowest = eivals(0,0).real();
    for (size_t i = 1, nRows = eivals.rows(); i < nRows; ++i) {
      double value = eivals(i,0).real();
      if(highest<value) highest = value;
      if(lowest>value) lowest = value;
    }
    return std::pair<double,double>(highest,lowest);
  }


std::vector<double> BCUT2D(const ROMol &m) {
  std::unique_ptr<ROMol> mol(MolOps::removeHs(m));
  // Atom elemno?  atomic mass? so many options...
  const std::string atom_prop = "atom_elemno";
  for(auto &atom: mol->atoms()) {
    atom->setProp<double>(atom_prop, atom->getMass());
  }
  // Gasteiger
  RDKit::computeGasteigerCharges(*mol, 12, true);
  // polarizability? - need model
  // slogp?  sasa?
  auto atom_bcut = BCUT2D(*mol, atom_prop);
  auto gasteiger = BCUT2D(*mol, common_properties::_GasteigerCharge);
  std::vector<double> res = {atom_bcut.first, atom_bcut.second, gasteiger.first, gasteiger.second};
  return res;
}
}
}

#endif
