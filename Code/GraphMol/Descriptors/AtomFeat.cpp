//
//  Copyright (C) 2020 Guillaume GODIN
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
#include <iostream>

#include "AtomFeat.h"

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <GraphMol/Atom.h>

#include <GraphMol/MolOps.h>
#include <cmath>
#include <vector>

namespace RDKit {

namespace Descriptors {

namespace {

std::vector<Atom::ChiralType> RS{Atom::CHI_TETRAHEDRAL_CW,
                                 Atom::CHI_TETRAHEDRAL_CCW, Atom::CHI_OTHER};
std::vector<std::string> Symbols{"B", "C",  "N",  "O", "S", "F", "Si",
                                 "P", "Cl", "Br", "I", "H", "*"};
std::vector<Atom::HybridizationType> HS{Atom::SP, Atom::SP2, Atom::SP3,
                                        Atom::SP3D, Atom::SP3D2};

void AtomFeatVector(const RDKit::Atom* atom, const ROMol* mol,
                    std::vector<double>& feats, bool addchiral) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(mol, "bad mol");

  if (addchiral) {
    feats.reserve(52);
  } else {
    feats.reserve(49);
  }

  // initiate ring info if not already done
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }

  // one hot atom symbols
  std::string s = atom->getSymbol();
  bool inlist = false;
  int indx = 0;
  for (auto ind : Symbols) {
    if (ind == s) {
      feats[indx] = 1;
      inlist = true;
    } else {
      feats[indx] = 0;
    }
    ++indx;
  }

  // write UNK type if not found in the symbol list
  feats[indx] = (inlist ? 0 : 1);
  ++indx;

  // one hot degree
  int d = atom->getDegree();
  for (int i = 0; i < 7; i++) {
    feats[indx] = d == i;
    ++indx;
  }

  Atom::HybridizationType hs = atom->getHybridization();
  // one hot hybridization type
  for (auto hsquery : HS) {
    feats[indx] = hs == hsquery;
    ++indx;
  }

  // one hot  Implicit Valence
  int IV = atom->getImplicitValence();
  for (int i = 0; i < 7; ++i) {
    feats[indx] = IV == i;
    ++indx;
  }

  // one hot  getFormalCharge
  int fc = atom->getFormalCharge();
  for (int i = -1; i < 2; i++) {
    feats[indx] = fc == i;
    ++indx;
  }

  // one hot ring size
  int atomid = atom->getIdx();
  for (unsigned int i = 3; i < 9; i++) {
    feats[indx] = mol->getRingInfo()->isAtomInRingOfSize(atomid, i);
    ++indx;
  }

  // Is aromatic
  feats[indx] = (atom->getIsAromatic());
  ++indx;

  // one hot  Total NumH
  unsigned int toth = atom->getTotalNumHs(false);
  for (unsigned int i = 0; i < 5; ++i) {
    feats[indx] = toth == i;
    ++indx;
  }

  // add numatoms
  feats[indx] = 1. / mol->getNumAtoms();
  ++indx;

  // put if here
  if (addchiral) {
    Atom::ChiralType rs = atom->getChiralTag();
    // one hot getChiralTag type
    for (auto rsquery : RS) {
      feats[indx] = rs == rsquery;
      ++indx;
    }
  }
}
}  // end of anonymous namespace

// entry point
void AtomFeatVect(const ROMol& mol, std::vector<double>& res, int atomid,
                  bool addchiral) {
  res.clear();
  if (addchiral) {
    res.resize(52);
  } else {
    res.resize(49);
  }

  AtomFeatVector(mol.getAtomWithIdx(atomid), &mol, res, addchiral);
}

}  // namespace Descriptors
}  // namespace RDKit