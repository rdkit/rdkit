//
//  Copyright (c) 2020, Brian Kelley
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Institue of Cancer Research.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
