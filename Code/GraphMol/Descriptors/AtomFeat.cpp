//
//  Copyright (c) 2020, Guillaume GODIN
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
//
// Adding ATOM FEATURES descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <iostream>

#include "AtomFeat.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Atom.h>

#include <GraphMol/MolOps.h>
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <deque>
#include <Eigen/Core>
#include <Eigen/QR>
#include <vector>


using namespace Eigen;

namespace RDKit {

namespace Descriptors {

namespace {

void AtomFeat1(const RDKit::Atom* atom, const ROMol* mol, std::vector <double> &feats,  bool addchiral) {

    // initiate ring info if not already done
    if( !mol->getRingInfo()->isInitialized() ) {
      RDKit::MolOps::findSSSR( *mol );
    }

    int atomid = atom->getIdx();
    int d = atom->getDegree();
    std::string s = atom->getSymbol();
    // one hot atom symbols
    bool inlist = false;
    int indx = 0;
    for (auto ind : Symbols) {
        if (ind == s) {
            feats[indx] = 1;
            inlist = true;
        } 
        else 
        {
            feats[indx] = 0;
        }
        indx+=1;
    }


    // write UNK type if not found in the symbol list
    feats[indx] = (inlist ? 0 : 1);
    indx+=1;

    // one hot degree
    for (int i=0;  i<7; i++) {
          feats[indx] = (d == i ? 1 : 0);
          indx+=1;
    }


    Atom::HybridizationType hs = atom->getHybridization();
    // one hot hybridization type
    for (auto hsquery:  HS) {
          feats[indx] = (hs == hsquery ? 1 : 0);
          indx+=1;
    }

    // one hot  Implicit Valence
    int IV = atom->getImplicitValence();
    for (int i=0;  i<7; i++) {
          feats[indx] = (IV == i ? 1 : 0);
          indx+=1;
    }

    // one hot  getFormalCharge
    int fc = atom->getFormalCharge();
    for (int i=-1;  i<2; i++) {
          feats[indx] = (fc == i ? 1 : 0);
          indx+=1;
    }

    // one hot ring size
    for (int i=3;  i<9; i++) {
          feats[indx] = (mol->getRingInfo()->isAtomInRingOfSize( atomid , i ) ? 1 : 0);
          indx+=1;
    }

    // Is aromatic
    feats[indx] = (atom->getIsAromatic());
    indx+=1;

    // one hot  Total NumH
    int toth = atom->getTotalNumHs(false);
    for (int i=0;  i<5; i++) {
          feats[indx] = (toth == i ? 1 : 0);
          indx+=1;
    }

    // add numatoms
    feats[indx] = 1./mol->getNumAtoms();
    indx+=1;


    // put if here
    if (addchiral) {
      Atom::ChiralType rs = atom->getChiralTag();
      // one hot getChiralTag type
      for (auto rsquery:  RS) {
          feats[indx] = (rs == rsquery ? 1 : 0);
          indx+=1;
      }
    }


}
}  // end of anonymous namespace

// entry point
void AtomFeat(const ROMol& mol, std::vector<double>& res, int atomid,  bool addchiral) {

  res.clear();
  if (addchiral) {
    res.resize(52);
  }
  else {
    res.resize(49);
  }

  AtomFeat1( mol.getAtomWithIdx(atomid),  &mol, res, addchiral);




}

}  // namespace Descriptors
}  // namespace RDKit