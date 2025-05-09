//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RWMol.h>
#include <cmath>
#include <GraphMol/MolOps.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include "EmbeddedFrag.h"
#include "DepictUtils.h"
#include "Templates.h"
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>
#include "RDDepictor.h"
#include <list>
#include <algorithm>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/Substruct/SubstructMatch.h>
constexpr double NEIGH_RADIUS = 2.5;

namespace RDDepict {
namespace {
// returns the atomic degree to be used for coordinate generation
unsigned int getDepictDegree(const RDKit::Atom *atom) {
  PRECONDITION(atom, "no atom");
  return atom->getDegree();
}
}  // end of anonymous namespace

EmbeddedFrag::EmbeddedFrag(unsigned int aid, const RDKit::ROMol *mol) {
  PRECONDITION(mol, "");
  PRECONDITION(aid < mol->getNumAtoms(), "");

  EmbeddedAtom eatm;
  eatm.aid = aid;
  RDGeom::Point2D org(0.0, 0.0);
  RDGeom::Point2D normal(1.0, 0.0);
  eatm.loc = org;
  eatm.normal = normal;
  eatm.angle = -1.0;
  eatm.ccw = true;
  eatm.neighs.clear();
  d_eatoms.clear();
  d_attachPts.clear();
  d_eatoms[aid] = eatm;
  d_done = false;
  dp_mol = mol;
  this->updateNewNeighs(aid);
}

EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
                           const RDKit::VECT_INT_VECT &fusedRings,
                           bool useRingTemplates) {
  PRECONDITION(mol, "");
  dp_mol = mol;
  d_eatoms.clear();
  d_attachPts.clear();
  this->embedFusedRings(fusedRings, useRingTemplates);
  d_done = false;
}

EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
                           const RDGeom::INT_POINT2D_MAP &coordMap) {
  // constructor of a case where the user specifies the coordinates for a
  // portion of the atoms in the molecule - we will use these coordinates
  // blindly without testing for any kind of correctness - user is GOD :)

  // we are not going to do much here simply add the atoms we have coordinates
  // for to this fragment; as a result this fragment may not be as ready to add
  // new neighbors etc. for the following reason.
  // - the user may have specified coords for only a part of the atoms in a
  //   fused ring systems
  // - once we use these coordinates we need to set up the atoms properly so
  //   that new neighbors can be added to them
  PRECONDITION(mol, "");
  dp_mol = mol;
  d_eatoms.clear();
  d_attachPts.clear();
  unsigned int na = mol->getNumAtoms();
  for (const auto &cri : coordMap) {
    unsigned int aid = cri.first;
    CHECK_INVARIANT(aid < na, "");
    EmbeddedAtom eatom(aid, cri.second);
    eatom.neighs.clear();
    eatom.df_fixed = true;
    d_eatoms[aid] = eatom;
    d_done = false;
  }
  this->setupNewNeighs();
  this->setupAttachmentPoints();
}

void EmbeddedFrag::computeNbrsAndAng(unsigned int aid,
                                     const RDKit::INT_VECT &doneNbrs) {
  //                                     const RDKit::ROMol *mol) {
  PRECONDITION(dp_mol, "");
  PRECONDITION(aid < dp_mol->getNumAtoms(), "");

  PRECONDITION(doneNbrs.size() >= 3, "");
  // we will find all the inter nbr angles, pick the one with the largest angle
  // make those neighbors the nbr1 and nbr2 of aid
  std::list<DOUBLE_INT_PAIR> anglePairs;
  double ang;
  for (auto nbi1 = doneNbrs.begin(); nbi1 != doneNbrs.end(); ++nbi1) {
    auto nbi3 = nbi1;
    for (auto nbi2 = nbi3++; nbi2 != doneNbrs.end(); ++nbi2) {
      ang = computeAngle(d_eatoms[aid].loc, d_eatoms[*nbi1].loc,
                         d_eatoms[*nbi2].loc);
      auto nbrPair = std::make_pair((*nbi1), (*nbi2));
      anglePairs.emplace_back(ang, nbrPair);
    }
  }
  anglePairs.sort([](auto pr1, auto pr2) { return pr1.first < pr2.first; });

  // more pain, more pain we unfortunately cannot right away pick the largest
  // angle - it is possible that we pick an angle that is in a fused ring - see
  // if I can explain this with a diagram
  //        _     _
  //       / B   C \                                this space
  //      /   \ /   \                               intentionally left blank
  //     |     A     |
  //     |     |     |
  //      \    D    /
  //       \_/   \_/
  //
  //  Let's say we are sitting on A with nbrs B, C, D - it is possible that we
  //  find ang(BAD) to be largest, but a new neighbor in this case will be added
  //  inside the ring We want to find ang(BAC) instead - which we will this do
  //  by checking that both our neighbors are not involved in more than one
  //  ring. Bridged systems - don't even go there
  auto winner = anglePairs.back();
  for (auto pr : boost::adaptors::reverse(anglePairs)) {
    if ((dp_mol->getRingInfo()->numAtomRings(pr.second.first) <= 1) &&
        (dp_mol->getRingInfo()->numAtomRings(pr.second.second) <= 1)) {
      winner = pr;
      break;
    }
  }

  auto winPair = winner.second;
  auto wnb1 = winPair.first;
  auto wnb2 = winPair.second;

  // now find the smallest angle that contains one of these nbrs
  int nb2 = -1, nb1 = -1;
  for (auto anglePair : anglePairs) {
    auto nbrPair = anglePair.second;
    if (wnb1 == nbrPair.first) {
      nb2 = wnb1;
      nb1 = nbrPair.second;
      break;
    } else if (wnb1 == nbrPair.second) {
      nb2 = wnb1;
      nb1 = nbrPair.first;
      break;
    } else if (wnb2 == nbrPair.first) {
      nb2 = wnb2;
      nb1 = nbrPair.second;
      break;
    } else if (wnb2 == nbrPair.second) {
      nb2 = wnb2;
      nb1 = nbrPair.first;
      break;
    }
  }

  // now find the rotation between nb1 and nb2
  auto wAng = winner.first;
  d_eatoms[aid].rotDir = rotationDir(d_eatoms[aid].loc, d_eatoms[nb1].loc,
                                     d_eatoms[nb2].loc, wAng);
  d_eatoms[aid].nbr1 = nb1;
  d_eatoms[aid].nbr2 = nb2;
  d_eatoms[aid].angle = 2 * M_PI - wAng;
}

// constructor to embed a cis/trans system
EmbeddedFrag::EmbeddedFrag(const RDKit::Bond *dblBond) {
  // Earlier embedding a cis/trans system meant to assign coordinates to the
  // atoms on the double bond as well as the neighboring atoms connected by the
  // single bond for which the cis/trans code has been specified. this causes
  // some ugliness in cases where these neighboring atoms are either part of a
  // different cis/trans system or a ring system. The function "merge" used to
  // deal with this ugliness. Now we will just embed the atoms on the double
  // bonds and mark at these atoms the direction in which the incoming single
  // bonds should go. Makes the merge function easier and address issue 171
  // simultaneously.
  PRECONDITION(dblBond, "");
  PRECONDITION(dblBond->getBondType() == RDKit::Bond::DOUBLE, "");
  auto stype = dblBond->getStereo();
  PRECONDITION(stype > RDKit::Bond::STEREOANY, "");
  const auto &nbrAtms = dblBond->getStereoAtoms();
  PRECONDITION(nbrAtms.size() == 2, "");
  dp_mol = &(dblBond->getOwningMol());

  auto begAtm = dblBond->getBeginAtomIdx();
  auto endAtm = dblBond->getEndAtomIdx();

  // the begin atom goes at the origin and the normal goes along -ve y-axis
  // to be rotate clock to add the cis/trans single bond
  EmbeddedAtom beatm;
  beatm.aid = begAtm;
  beatm.loc = RDGeom::Point2D(0.0, 0.0);
  beatm.nbr1 = endAtm;

  beatm.normal = RDGeom::Point2D(0.0, -1.0);
  beatm.ccw = false;
  beatm.CisTransNbr = nbrAtms[0];
  d_eatoms[begAtm] = beatm;

  // the end atom goes on the x-axis
  EmbeddedAtom eeatm;
  eeatm.aid = endAtm;
  eeatm.loc = RDGeom::Point2D(BOND_LEN, 0.0);
  eeatm.nbr1 = begAtm;
  eeatm.CisTransNbr = nbrAtms[1];
  if (stype == RDKit::Bond::STEREOZ || stype == RDKit::Bond::STEREOCIS) {
    eeatm.normal = RDGeom::Point2D(0.0, -1.0);
    eeatm.ccw = true;
  } else {
    eeatm.normal = RDGeom::Point2D(0.0, 1.0);
    eeatm.ccw = false;
  }
  d_eatoms[endAtm] = eeatm;

  d_done = false;
}

int EmbeddedFrag::findNumNeigh(const RDGeom::Point2D &pt, double radius) {
  // find the number of atoms in the current embedded system that are within
  // 'radius' of the specified point
  int res = 0;
  for (const auto &efi : d_eatoms) {
    const auto &rloc = efi.second.loc;
    if ((rloc - pt).length() < radius) {
      ++res;
    }
  }
  return res;
}

void EmbeddedFrag::updateNewNeighs(
    unsigned int aid) {  //, const RDKit::ROMol *mol) {
  PRECONDITION(dp_mol, "");

  d_eatoms[aid].neighs.clear();
  RDKit::INT_VECT hIndices;
  for (const auto nbr : dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(aid))) {
    if (d_eatoms.find(nbr->getIdx()) == d_eatoms.end()) {
      if (dp_mol->getAtomWithIdx(nbr->getIdx())->getAtomicNum() != 1) {
        d_eatoms[aid].neighs.push_back(nbr->getIdx());
      } else {
        hIndices.push_back(nbr->getIdx());
      }
    }
  }
  d_eatoms[aid].neighs.insert(d_eatoms[aid].neighs.end(), hIndices.begin(),
                              hIndices.end());

  auto deg = getDepictDegree(dp_mol->getAtomWithIdx(aid));
  // order the neighbors by their CIPranks, if the number is between > 0 but
  // less than 3
  if ((d_eatoms[aid].neighs.size() > 0) &&
      ((deg < 4) || (d_eatoms[aid].neighs.size() < 3))) {
    d_eatoms[aid].neighs = rankAtomsByRank(*dp_mol, d_eatoms[aid].neighs);
  } else if ((deg >= 4) && (d_eatoms[aid].neighs.size() >= 3)) {
    // now if we have more more than 2 neighbors change the order so that atoms
    // with the highest rank fall on opposite sides of each other
    d_eatoms[aid].neighs = setNbrOrder(aid, d_eatoms[aid].neighs, *dp_mol);
  }

  if (d_eatoms[aid].neighs.size() > 0) {
    if (std::find(d_attachPts.begin(), d_attachPts.end(),
                  static_cast<int>(aid)) == d_attachPts.end()) {
      d_attachPts.push_back(aid);
    }
  }
}

void EmbeddedFrag::setupNewNeighs() {  // const RDKit::ROMol *mol) {
  PRECONDITION(dp_mol, "");

  d_attachPts.clear();
  for (const auto &eci : d_eatoms) {
    this->updateNewNeighs(eci.first);
  }
  // arrange the d_attachPts so that they are traversed in the order of CIPRanks
  d_attachPts = rankAtomsByRank(*dp_mol, d_attachPts);
}

int EmbeddedFrag::findNeighbor(
    unsigned int aid) {  //, const RDKit::ROMol *mol) {
  PRECONDITION(dp_mol, "");

  for (const auto nbr : dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(aid))) {
    if (d_eatoms.find(nbr->getIdx()) != d_eatoms.end()) {
      return nbr->getIdx();
    }
  }
  return -1;
}

void EmbeddedFrag::setupAttachmentPoints() {
  // now for points that new atoms will be added to later on we need to do some
  // setup
  for (auto dai : d_attachPts) {
    // find the neighbors that are already embedded for each of these atoms
    RDKit::INT_VECT doneNbrs;
    const auto &enbrs = d_eatoms[dai].neighs;
    for (const auto nbrAtom :
         dp_mol->atomNeighbors(dp_mol->getAtomWithIdx(dai))) {
      if (std::find(enbrs.begin(), enbrs.end(),
                    static_cast<int>(nbrAtom->getIdx())) == enbrs.end()) {
        // we found a neighbor that is part of this embedded system
        doneNbrs.push_back(nbrAtom->getIdx());
      }
    }
    if (doneNbrs.empty()) {
      d_eatoms[dai].normal = RDGeom::Point2D(1., 0.);
      d_eatoms[dai].angle = -1.;
    } else if (doneNbrs.size() == 1) {
      auto nbid = doneNbrs.front();
      d_eatoms[dai].nbr1 = nbid;
      d_eatoms[dai].normal =
          computeNormal(d_eatoms[dai].loc, d_eatoms[nbid].loc);
    } else if (doneNbrs.size() == 2) {
      auto nb1 = doneNbrs[0];
      auto nb2 = doneNbrs[1];
      d_eatoms[dai].nbr1 = nb1;
      d_eatoms[dai].nbr2 = nb2;
      d_eatoms[dai].angle =
          computeAngle(d_eatoms[dai].loc, d_eatoms[nb1].loc, d_eatoms[nb2].loc);
    } else if (doneNbrs.size() >= 3) {
      // this is a pain - delegate it to a utility function
      this->computeNbrsAndAng(dai, doneNbrs);
    }
  }
}

// check if the stereochemistry of the template matches the stereochemistry of
// the molecule
static bool checkStereoChemistry(const RDKit::ROMol &mol,
                                 const RDKit::ROMol &template_mol,
                                 RDKit::MatchVectType match) {
  for (auto bond : mol.bonds()) {
    if (bond->getBondType() != RDKit::Bond::DOUBLE ||
        bond->getStereo() == RDKit::Bond::STEREOANY ||
        bond->getStereo() == RDKit::Bond::STEREONONE) {
      continue;
    }
    // get the four atoms around the double bond
    auto neighbors = bond->getStereoAtoms();
    if (neighbors.size() != 2) {
      continue;
    }
    int atom1_neighbor1 = neighbors[0];
    int atom2_neighbor1 = neighbors[1];
    int atom1 = bond->getBeginAtomIdx();
    int atom2 = bond->getEndAtomIdx();

    // now get the other two atoms that are not part of the double bond (if any)
    int atom1_neighbor2 = -1;
    int atom2_neighbor2 = -1;
    if (mol.getAtomWithIdx(atom1)->getDegree() > 2) {
      for (auto neighbor : mol.atomNeighbors(mol.getAtomWithIdx(atom1))) {
        if (static_cast<int>(neighbor->getIdx()) != atom1_neighbor1 &&
            static_cast<int>(neighbor->getIdx()) != atom2) {
          atom1_neighbor2 = neighbor->getIdx();
          break;
        }
      }
    }
    if (mol.getAtomWithIdx(atom2)->getDegree() > 2) {
      for (auto neighbor : mol.atomNeighbors(mol.getAtomWithIdx(atom2))) {
        if (static_cast<int>(neighbor->getIdx()) != atom2_neighbor1 &&
            static_cast<int>(neighbor->getIdx()) != atom1) {
          atom2_neighbor2 = neighbor->getIdx();
          break;
        }
      }
    }

    // find the template atoms that correspond to the four atoms
    int template_atom1 = -1;
    int template_atom2 = -1;
    int template_atom1_neighbor1 = -1;
    int template_atom1_neighbor2 = -1;
    int template_atom2_neighbor1 = -1;
    int template_atom2_neighbor2 = -1;
    for (auto &[template_aidx, rs_aidx] : match) {
      if (rs_aidx == atom1) {
        template_atom1 = template_aidx;
      } else if (rs_aidx == atom2) {
        template_atom2 = template_aidx;
      } else if (rs_aidx == atom1_neighbor1) {
        template_atom1_neighbor1 = template_aidx;
      } else if (rs_aidx == atom2_neighbor1) {
        template_atom2_neighbor1 = template_aidx;
      } else if (rs_aidx == atom1_neighbor2) {
        template_atom1_neighbor2 = template_aidx;
      } else if (rs_aidx == atom2_neighbor2) {
        template_atom2_neighbor2 = template_aidx;
      }
    }

    // there's a chance that the atoms controlling the double bond stereochem in
    // the molecule are not the atoms that matched to the template, handle that
    // here by swapping to the other atom
    bool swapStereo = false;
    if (template_atom1_neighbor1 == -1) {
      template_atom1_neighbor1 = template_atom1_neighbor2;
      swapStereo = !swapStereo;
    }
    if (template_atom2_neighbor1 == -1) {
      template_atom2_neighbor1 = template_atom2_neighbor2;
      swapStereo = !swapStereo;
    }

    if (template_atom1 == -1 || template_atom2 == -1 ||
        template_atom1_neighbor1 == -1 || template_atom2_neighbor1 == -1) {
      return false;
    }

    const auto &conf = template_mol.getConformer();
    const auto &atom1_loc = conf.getAtomPos(template_atom1);
    const auto &atom2_loc = conf.getAtomPos(template_atom2);
    const auto &atom1_neighbor_loc = conf.getAtomPos(template_atom1_neighbor1);
    const auto &atom2_neighbor_loc = conf.getAtomPos(template_atom2_neighbor1);
    // check if the two neighbors are on the same side of the bond
    const auto v12 = atom1_neighbor_loc - atom1_loc;
    const auto v42 = atom2_neighbor_loc - atom1_loc;
    const auto v32 = atom2_loc - atom1_loc;
    auto cross1 = v32.x * v12.y - v32.y * v12.x;
    auto cross2 = v32.x * v42.y - v32.y * v42.x;
    bool is_cis = cross1 * cross2 > 0;
    if (swapStereo) {
      is_cis = !is_cis;
    }
    if (is_cis != (bond->getStereo() == RDKit::Bond::STEREOZ ||
                   bond->getStereo() == RDKit::Bond::STEREOCIS)) {
      return false;
    }
  }
  return true;
}

bool EmbeddedFrag::matchToTemplate(const RDKit::INT_VECT &ringSystemAtoms,
                                   unsigned int ring_count) {
  CoordinateTemplates &coordinate_templates =
      CoordinateTemplates::getRingSystemTemplates();

  // only look for an exact match to the ring system because our method of
  // completing rings from a template isn't reliably better than not using
  // a template at all
  if (!coordinate_templates.hasTemplateOfSize(ringSystemAtoms.size())) {
    return false;
  }

  // make a mol out of the induced subgraph using the ring system atoms
  RDKit::RWMol rs_mol(*dp_mol, true);

  boost::dynamic_bitset<> rs_atoms(dp_mol->getNumAtoms());
  for (auto aidx : ringSystemAtoms) {
    rs_atoms.set(aidx);
  }

  constexpr int DUMMY_ATOMIC_NUM = 200;
  for (auto &at : rs_mol.atoms()) {
    if (!rs_atoms.test(at->getIdx())) {
      at->setAtomicNum(DUMMY_ATOMIC_NUM);
    }
  }
  auto numBonds = rs_mol.getNumBonds();
  for (auto bnd : rs_mol.bonds()) {
    if (!rs_atoms.test(bnd->getBeginAtomIdx()) ||
        !rs_atoms.test(bnd->getEndAtomIdx())) {
      --numBonds;
    }
  }

  // find template that this mol matches to, if any
  RDKit::MatchVectType match;
  std::shared_ptr<RDKit::ROMol> template_mol(nullptr);
  for (const auto &mol :
       coordinate_templates.getMatchingTemplates(ringSystemAtoms.size())) {
    // To reduce how often we have to do substructure matches, check ring info
    // and bond count first
    if (mol->getNumBonds() != numBonds) {
      continue;
    } else if (mol->getRingInfo()->numRings() != ring_count) {
      continue;
    }
    // also check if the mol atoms have the same connectivity as the template
#ifdef _MSC_VER
    // MSVC++ doesn't like implicitly capturing constexpr variables, this is a
    // bug
    auto degreeCounts = [DUMMY_ATOMIC_NUM](const RDKit::ROMol &mol) {
#else
    // clang generates warnings if you explicitly capture a constexpr variable
    auto degreeCounts = [](const RDKit::ROMol &mol) {
#endif
      std::array<int, 5> degrees_count({0, 0, 0, 0, 0});
      for (auto atom : mol.atoms()) {
        if (atom->getAtomicNum() == DUMMY_ATOMIC_NUM) {
          continue;
        }
        auto degree = 0u;
        for (auto nbr : mol.atomNeighbors(atom)) {
          if (nbr->getAtomicNum() != DUMMY_ATOMIC_NUM) {
            ++degree;
            if (degree == 4) {
              break;
            }
          }
        }
        degrees_count[degree]++;
      }
      return degrees_count;
    };
    if (degreeCounts(rs_mol) != degreeCounts(*mol)) {
      continue;
    }
    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1;
    auto matches = RDKit::SubstructMatch(rs_mol, *mol, params);
    if (!matches.empty()) {
      if (checkStereoChemistry(rs_mol, *mol, matches[0])) {
        match = matches[0];
        template_mol = mol;
        break;
      }
    }
  }
  if (!template_mol) {
    return false;
  }

  // copy over new coordinates
  const auto &conf = template_mol->getConformer();
  for (auto &[template_aidx, rs_aidx] : match) {
    EmbeddedAtom new_at(rs_aidx, conf.getAtomPos(template_aidx));
    new_at.df_fixed = true;
    d_eatoms.emplace(rs_aidx, new_at);
  }
  this->setupNewNeighs();
  this->setupAttachmentPoints();
  return true;
}

// find any atoms in the ring that are in trans double bonds
// and mirror them into the ring
static void mirrorTransRingAtoms(const RDKit::ROMol &mol,
                                 const RDKit::INT_VECT &ring,
                                 RDGeom::INT_POINT2D_MAP &coords) {
  // a nice place for C++23 generator coroutines...
  RDKit::INT_VECT transRingAtoms;
  for (size_t i = 0; i < ring.size(); ++i) {
    const auto atom1 = ring[i];
    const auto atom2 = ring[(i + 1) % ring.size()];
    const auto bond = mol.getBondBetweenAtoms(atom1, atom2);
    if (bond->getBondType() != RDKit::Bond::DOUBLE) {
      continue;
    }
    const auto stype = bond->getStereo();
    if (stype <= RDKit::Bond::STEREOANY) {
      continue;
    }

    // We care about bonds that are trans with respect to this ring
    const auto &neighbors = bond->getStereoAtoms();
    if (neighbors.size() != 2) {
      continue;
    }
    const auto leftIsIn =
        std::find(ring.begin(), ring.end(), neighbors[0]) != ring.end();
    const auto rightIsIn =
        std::find(ring.begin(), ring.end(), neighbors[1]) != ring.end();
    bool isTrans = false;
    if (stype == RDKit::Bond::STEREOTRANS || stype == RDKit::Bond::STEREOE) {
      if (leftIsIn == rightIsIn) {
        // trans, both neighbors in the ring (or both out)
        isTrans = true;
      }
    } else if (leftIsIn != rightIsIn) {
      // cis, but one of the neighbors is outside the ring
      isTrans = true;
    }
    if (!isTrans) {
      continue;
    }

    // Mirror one atom in each trans bond across the line defined by its two
    // neighbors. This bumps it into the ring
    const auto left = ring[(i + ring.size() - 1) % ring.size()];
    const auto right = atom2;

    const auto last = coords[left];
    const auto ref = coords[right];
    const auto interest = coords[atom1];
    const auto d = last - ref;
    const double a = (d.x * d.x - d.y * d.y) / d.dotProduct(d);
    const double b = 2 * d.x * d.y / d.dotProduct(d);
    const double x =
        a * (interest.x - ref.x) + b * (interest.y - ref.y) + ref.x;
    const double y =
        b * (interest.x - ref.x) - a * (interest.y - ref.y) + ref.y;
    coords[atom1] = RDGeom::Point2D(x, y);
  }
}

//
// NOTE: the individual rings in fusedRings must appear in traversal order.
//    This is what is provided by the current ring-finding code.
//
void EmbeddedFrag::embedFusedRings(const RDKit::VECT_INT_VECT &fusedRings,
                                   bool useRingTemplates) {
  PRECONDITION(dp_mol, "");
  // Look for a template for the whole system. Failing that simplify the system
  // to a set of core atoms and  look for a template for those. If that fails,
  // start from a single ring. Then add rings one by one

  RDKit::INT_VECT funion;
  // look for a template that matches the entire fused ring system
  if (useRingTemplates && fusedRings.size() > 1) {
    RDKit::Union(fusedRings, funion);
    bool found_template = matchToTemplate(funion, fusedRings.size());
    if (found_template) {
      // we are done
      return;
    }
  }
  std::vector<RDGeom::INT_POINT2D_MAP> coords;
  coords.reserve(fusedRings.size());

  for (const auto &ring : fusedRings) {
    auto ring_coords = embedRing(ring);
    mirrorTransRingAtoms(*dp_mol, ring, ring_coords);
    coords.push_back(ring_coords);
  }
  RDKit::INT_VECT doneRings;

  if (useRingTemplates) {
    RDKit::INT_VECT coreRingsIds;
    auto coreRings = findCoreRings(fusedRings, coreRingsIds, *dp_mol);
    if (coreRings.size() > 1 && coreRings.size() < fusedRings.size()) {
      // look for a template that matches the core ring system
      RDKit::Union(coreRings, funion);
      bool found_template = matchToTemplate(funion, coreRings.size());
      if (found_template) {
        doneRings = coreRingsIds;
      }
    }
  }

  // if not embed find a ring as a starting point
  if (doneRings.empty()) {
    // FIX for issue 197
    // find the ring with the max substituents
    // If there are multiple pick the largest
    auto firstRingId = pickFirstRingToEmbed(*dp_mol, fusedRings);

    this->initFromRingCoords(fusedRings[firstRingId], coords[firstRingId]);
    doneRings.push_back(firstRingId);
  }
  RDKit::Union(fusedRings, funion);
  // now loop over the remaining rings and attach them one at a time
  // the order is determined by how many atoms a ring has in common with
  // the atoms already embedded
  while (d_eatoms.size() < funion.size()) {  // ) {
    int nextId;
    // we will take the ring with maximum number of common atoms with
    // with atoms already done
    auto commonAtomIds = findNextRingToEmbed(doneRings, fusedRings, nextId);

    RDGeom::Transform2D trans;
    EmbeddedFrag embRing;
    embRing.initFromRingCoords(fusedRings[nextId], coords[nextId]);
    RDKit::INT_VECT pinAtoms;
    // REVIEW: using the average position of the shared atoms and the
    // centroid vector, we can make this a single case.
    if (commonAtomIds.size() == 1) {
      trans.assign(this->computeOneAtomTrans(commonAtomIds[0], embRing));
      embRing.Transform(trans);
      pinAtoms.push_back(commonAtomIds.front());
    } else {
      // if the common atoms form a chain they are going to be in order - we try
      // to do that in findNextRingToEmbed we will therefore try to use the last
      // and the first atoms in the chain to fuse the rings - will hopefully fix
      // issue 177
      auto aid1 = commonAtomIds.front();
      auto aid2 = commonAtomIds.back();
      pinAtoms.push_back(aid1);
      pinAtoms.push_back(aid2);
      trans.assign(this->computeTwoAtomTrans(aid1, aid2, coords[nextId]));
      embRing.Transform(trans);
      reflectIfNecessaryDensity(embRing, aid1, aid2);
    }
    this->mergeRing(embRing, commonAtomIds.size(), pinAtoms);
    doneRings.push_back(nextId);
  }
}

RDGeom::Transform2D EmbeddedFrag::computeOneAtomTrans(
    unsigned int commAid, const EmbeddedFrag &other) {
  // find the coordinates for the same atom in the embedded system
  auto rcr = d_eatoms[commAid].loc;

  // find the coordinate for the same atom in the other system
  const auto &oeatm = other.GetEmbeddedAtom(commAid);
  auto ccr = oeatm.loc;
  auto onb1 = oeatm.nbr1;
  auto onb2 = oeatm.nbr2;
  CHECK_INVARIANT((onb1 >= 0) && (onb2 >= 0), "");
  auto midPt = other.GetEmbeddedAtom(onb1).loc;
  midPt += other.GetEmbeddedAtom(onb2).loc;
  midPt *= 0.5;

  // get the coordinates for the neighboring atoms
  auto nb1 = d_eatoms[commAid].nbr1;
  auto nb2 = d_eatoms[commAid].nbr2;
  auto nbp1 = d_eatoms[nb1].loc;
  auto nbp2 = d_eatoms[nb2].loc;

  auto ang = d_eatoms[commAid].angle;
  auto largestAngle = 2 * M_PI - ang;

  auto bpt = computeBisectPoint(rcr, largestAngle, nbp1, nbp2);

  // now that we have the bisect point compute the transform that will take ccr
  // to coincide with rcr and the mid point between the neighbors of ccr to fall
  // on the line from rcr to bpt
  RDGeom::Transform2D trans;
  trans.SetTransform(rcr, bpt, ccr, midPt);
  return trans;
}

RDGeom::Transform2D EmbeddedFrag::computeTwoAtomTrans(
    unsigned int aid1, unsigned int aid2,
    const RDGeom::INT_POINT2D_MAP &nringCor) {
  CHECK_INVARIANT(d_eatoms.find(aid1) != d_eatoms.end(), "");
  CHECK_INVARIANT(d_eatoms.find(aid2) != d_eatoms.end(), "");

  // this is an easier thing to do than computeOneAtomTrans
  // we know that there are at least two atoms in common between the new ring
  // and the rings that have already been embedded.
  //
  // we are going to simply use the first two atoms on the commIds list and
  // use those to compute a transforms
  const auto &loc1 = nringCor.at(aid1);
  const auto &loc2 = nringCor.at(aid2);

  // get the coordinates for the same atoms in the already embedded ring system
  const auto &ref1 = d_eatoms.at(aid1).loc;
  const auto &ref2 = d_eatoms.at(aid2).loc;
  RDGeom::Transform2D trans;
  trans.SetTransform(ref1, ref2, loc1, loc2);
  return trans;
}

void EmbeddedFrag::Reflect(const RDGeom::Point2D &loc1,
                           const RDGeom::Point2D &loc2) {
  for (auto &ei : d_eatoms) {
    ei.second.Reflect(loc1, loc2);
  }
}

void EmbeddedFrag::reflectIfNecessaryCisTrans(EmbeddedFrag &embFrag,
                                              unsigned int ctCase,
                                              unsigned int aid1,
                                              unsigned int aid2) {
  // ok this is a cis/trans case - we may have violated the cis/trans
  // specification
  // so lets try to correct it with a reflection
  const auto &p1Loc = d_eatoms[aid1].loc;
  RDGeom::Point2D rAtmLoc, p1norm;
  if (ctCase == 1) {
    // embObj is the cis/trans case - find the normal at aid1 - this should tell
    // us where the ring single bond in the cis/trans system should have gone
    p1norm = embFrag.d_eatoms[aid1].normal;
    auto ringAtm = embFrag.d_eatoms[aid1].CisTransNbr;
    if (d_eatoms.find(ringAtm) != d_eatoms.end()) {
      rAtmLoc = d_eatoms[ringAtm].loc;
    } else {
      // FIX: this is a work-around arising from issue 3135833
      BOOST_LOG(rdWarningLog) << "Warning: stereochemistry around double bond "
                                 "may be incorrect in depiction."
                              << std::endl;
      return;
    }
  } else {
    // this is the cis/trans object
    p1norm = d_eatoms[aid1].normal;
    auto ringAtm = d_eatoms[aid1].CisTransNbr;
    rAtmLoc = embFrag.d_eatoms[ringAtm].loc;
  }
  rAtmLoc -= p1Loc;
  auto dot = rAtmLoc.dotProduct(p1norm);
  auto p2Loc = d_eatoms[aid2].loc;
  if (dot < 0.0) {
    embFrag.Reflect(p1Loc, p2Loc);
  }
}

void EmbeddedFrag::reflectIfNecessaryThirdPt(EmbeddedFrag &embFrag,
                                             unsigned int aid1,
                                             unsigned int aid2,
                                             unsigned int aid3) {
  const auto &pt1 = d_eatoms[aid1].loc;
  const auto &pt2 = d_eatoms[aid2].loc;

  auto normal = pt2;
  normal -= pt1;
  normal.rotate90();

  const auto oth3 = embFrag.GetEmbeddedAtom(aid3).loc - pt1;
  const auto pt3 = d_eatoms[aid3].loc - pt1;

  auto dot1 = normal.dotProduct(pt3);
  auto dot2 = normal.dotProduct(oth3);
  if (dot1 * dot2 < 0.0) {
    // the third atom is on either sides of the line between aid1 and aid2 in
    // the two fragment - let us reflect to correct it
    embFrag.Reflect(pt1, pt2);
  }
}

void EmbeddedFrag::reflectIfNecessaryDensity(EmbeddedFrag &embFrag,
                                             unsigned int aid1,
                                             unsigned int aid2) {
  // ok we will do this the new way by measuring a density function
  const auto &pin1 = d_eatoms[aid1].loc;
  const auto &pin2 = d_eatoms[aid2].loc;
  double densityNormal = 0.0;
  double densityReflect = 0.0;
  for (const auto &oci : embFrag.GetEmbeddedAtoms()) {
    if (d_eatoms.find(oci.first) == d_eatoms.end()) {
      auto loc1 = oci.second.loc;
      auto rloc1 = reflectPoint(loc1, pin1, pin2);
      for (const auto &tci : d_eatoms) {
        auto t1 = tci.second.loc;
        t1 -= loc1;
        auto td = t1.length();
        auto rt1 = tci.second.loc;
        rt1 -= rloc1;
        auto rtd = rt1.length();
        if (td > 1.0e-3) {
          densityNormal += (1.0 / td);
        } else {
          densityNormal += 1000.0;
        }
        if (rtd > 1.0e-3) {
          densityReflect += (1.0 / rtd);
        } else {
          densityReflect += 1000.0;
        }
      }
    }
  }
  if (densityNormal - densityReflect > 1.0e-4) {
    embFrag.Reflect(pin1, pin2);
  }
}

void EmbeddedFrag::initFromRingCoords(const RDKit::INT_VECT &ring,
                                      const RDGeom::INT_POINT2D_MAP &nringMap) {
  double largestAngle = M_PI * (1 - (2.0 / ring.size()));
  auto prev = ring.back();
  unsigned int cnt = 0;
  for (auto ai : ring) {
    EmbeddedAtom eatm;
    eatm.loc = nringMap.at(ai);
    eatm.aid = ai;
    eatm.angle = largestAngle;
    eatm.nbr1 = prev;
    if (cnt) {
      d_eatoms[prev].nbr2 = ai;
    }
    d_eatoms[ai] = eatm;
    prev = ai;
    cnt++;
  }
  d_eatoms[prev].nbr2 = ring.front();
}

void EmbeddedFrag::mergeRing(const EmbeddedFrag &embRing, unsigned int nCommon,
                             const RDKit::INT_VECT &pinAtoms) {
  const auto &oatoms = embRing.GetEmbeddedAtoms();
  for (const auto &ori : oatoms) {
    auto aid = ori.first;
    if (d_eatoms.find(aid) == d_eatoms.end()) {
      d_eatoms[aid] = ori.second;
    } else {
      // update the neighbor only on atoms that were used to compute the
      // transform to merge the and only if the two are the only common atoms
      // i.e. we are doing bridged systems we will leave the nbrs untouched
      if (nCommon <= 2) {
        if (std::find(pinAtoms.begin(), pinAtoms.end(), aid) !=
            pinAtoms.end()) {
          d_eatoms[aid].angle += ori.second.angle;
          if (d_eatoms[aid].nbr1 == ori.second.nbr1) {
            d_eatoms[aid].nbr1 = ori.second.nbr2;
          } else if (d_eatoms[aid].nbr1 == ori.second.nbr2) {
            d_eatoms[aid].nbr1 = ori.second.nbr1;
          } else if (d_eatoms[aid].nbr2 == ori.second.nbr1) {
            d_eatoms[aid].nbr2 = ori.second.nbr2;
          } else if (d_eatoms[aid].nbr2 == ori.second.nbr2) {
            d_eatoms[aid].nbr2 = ori.second.nbr1;
          }
        }
      }
    }
  }
}

void EmbeddedFrag::addNonRingAtom(unsigned int aid, unsigned int toAid) {
  // const RDKit::ROMol *mol) {
  PRECONDITION(dp_mol, "");
  // check that aid does not belong the embedded fragment yet
  PRECONDITION(d_eatoms.find(aid) == d_eatoms.end(), "");
  // and that toAid is already in the embedded system
  PRECONDITION(d_eatoms.find(toAid) != d_eatoms.end(), "");
  if (d_eatoms[toAid].angle > 0.0) {
    addAtomToAtomWithAng(aid, toAid);
  } else {
    addAtomToAtomWithNoAng(aid, toAid);
  }
  // remove aid from the neighbor list of toAid
  d_eatoms[toAid].neighs.erase(std::remove(d_eatoms[toAid].neighs.begin(),
                                           d_eatoms[toAid].neighs.end(),
                                           static_cast<int>(aid)));
  this->updateNewNeighs(aid);
}

void EmbeddedFrag::addAtomToAtomWithAng(unsigned int aid, unsigned int toAid) {
  const auto &refAtom = d_eatoms[toAid];
  auto refLoc = refAtom.loc;
  RDGeom::Point2D origin(0.0, 0.0);
  PRECONDITION(refAtom.angle > 0.0, "");

  // we are adding to either to a ring atom or an atom to which we added at
  // least one substituent previously

  // determine the angle at which we want to add the new atom based on the
  // number of remaining substituents
  auto nnbr = refAtom.neighs.size();
  double remAngle = 2 * M_PI - refAtom.angle;
  auto currAngle = remAngle / (1 + nnbr);
  d_eatoms[toAid].angle += currAngle;

  const auto &nb1 = d_eatoms.at(refAtom.nbr1).loc;
  const auto &nb2 = d_eatoms.at(refAtom.nbr2).loc;
  if (d_eatoms[toAid].rotDir == 0) {
    d_eatoms[toAid].rotDir = rotationDir(refLoc, nb1, nb2, remAngle);
  }

  currAngle *= d_eatoms[toAid].rotDir;

  RDGeom::Transform2D rtrans;
  rtrans.SetTransform(refLoc, currAngle);
  auto currLoc = nb2;
  rtrans.TransformPoint(currLoc);
  if (fabs(remAngle) - M_PI < 1e-3) {
    auto currLoc2 = nb2;
    rtrans.SetTransform(refLoc, -currAngle);
    rtrans.TransformPoint(currLoc2);
    if (findNumNeigh(currLoc, 0.5) > findNumNeigh(currLoc2, 0.5)) {
      currLoc = currLoc2;
      currAngle *= -1;
    } else {
      rtrans.SetTransform(refLoc, currAngle);
    }
  }

  // set the neighbors for the current point
  d_eatoms[toAid].nbr2 = aid;

  EmbeddedAtom eatm;
  eatm.aid = aid;
  eatm.loc = currLoc;
  eatm.nbr1 = toAid;
  eatm.angle = -1.0;

  // now compute the normal at this atom - which gives the direction in which we
  // want to add the next atom. We will go in the direction that seem to be
  // least explored
  auto tpt = currLoc - refLoc;
  RDGeom::Point2D norm(-tpt.y, tpt.x);
  auto tp1 = currLoc + norm;
  auto tp2 = currLoc - norm;

  auto nccw = findNumNeigh(
      tp1, NEIGH_RADIUS);  // number of neighbors if we go counter-clockwise
  auto ncw = findNumNeigh(
      tp2, NEIGH_RADIUS);  // number of neighbors if we go clockwise

  norm.normalize();
  if (nccw < ncw) {
    eatm.normal = norm;
    eatm.ccw = false;
  } else {
    eatm.normal = (-norm);
    eatm.ccw = true;
  }

  d_eatoms[aid] = eatm;
}

void EmbeddedFrag::addAtomToAtomWithNoAng(unsigned int aid,
                                          unsigned int toAid) {
  PRECONDITION(dp_mol, "");
  const auto &refAtom = d_eatoms.at(toAid);
  PRECONDITION(refAtom.angle <= 0.0, "");
  const auto &refLoc = refAtom.loc;
  RDGeom::Point2D origin(0.0, 0.0);
  auto refAtomCCW = refAtom.ccw;

  // -----------------------------------------------------------------------
  // we are adding to a non-ring atom,
  // the direction in which we add the new atom matters here
  auto currLoc = refAtom.normal;
  if (refAtom.CisTransNbr >= 0) {
    // ok this atom is part of a cis/trans dbl bond
    if (static_cast<unsigned int>(refAtom.CisTransNbr) != aid) {
      // but we are note adding the single bond atom to which the cis/trans
      // specification was made, in this case reverse the normal and the ccw
      refAtomCCW = !refAtomCCW;
      currLoc *= -1.0;
    }
  }

  CHECK_INVARIANT(currLoc.lengthSq() > 1.0e-8, "");

  // find out what angle we want to add bond at
  const auto atm = dp_mol->getAtomWithIdx(toAid);
  auto deg = getDepictDegree(atm);

  auto angle = computeSubAngle(deg, atm->getHybridization());

  // update the current atom we already have a nbr1 set on the current atom
  // update the angle etc d_eatoms[toAid].nbr2 = aid;
  bool flipNorm = false;
  if (d_eatoms[toAid].nbr1 >= 0) {
    d_eatoms[toAid].angle = angle;

    d_eatoms[toAid].nbr2 = aid;
  } else {
    // ------------------
    // We'll be here for the first atom in a system with no rings, we have
    // nothing
    // else set up, so we will deal with this case carefully.
    //  - if the angle is 120 deg we will add the first atom at 30 deg angle to
    //  the x-axis
    //  - for any other angle we will use the x-axis to add the new atom
    //  - we will set the normal perpendicular to this first bond in the counter
    //  clockwise direction
    //
    // RDGeom::Point2D norm;

    auto norm = d_eatoms.at(toAid).normal;
    RDGeom::Transform2D rtrans;
    rtrans.SetTransform(origin, angle);
    rtrans.TransformPoint(norm);
    d_eatoms[toAid].normal = norm;
    d_eatoms[toAid].nbr1 = aid;
    flipNorm = true;
  }

  angle -= M_PI / 2;
  if (!refAtomCCW) {
    // we want to rotate clockwise
    angle *= -1.0;
  }

  RDGeom::Transform2D trans;
  trans.SetTransform(origin, angle);
  trans.TransformPoint(currLoc);
  currLoc *= BOND_LEN;
  currLoc += refLoc;

  // now compute the normal at this new point for the next addition
  auto tpt = refLoc - currLoc;
  // This is the lazy man's rotation by 90 degrees about the origin:
  RDGeom::Point2D norm(-tpt.y, tpt.x);
  if (refAtomCCW ^ flipNorm) {
    norm *= -1.0;
  }
  norm.normalize();
  EmbeddedAtom eatm;
  eatm.loc = currLoc;
  eatm.normal = norm;
  eatm.nbr1 = toAid;

  eatm.angle = -1.0;

  eatm.ccw = (!refAtomCCW) ^ flipNorm;
  d_eatoms[aid] = eatm;
}

RDKit::INT_VECT EmbeddedFrag::findCommonAtoms(const EmbeddedFrag &efrag2) {
  RDKit::INT_VECT res;
  for (auto eri1 : this->GetEmbeddedAtoms()) {
    for (auto eri2 : efrag2.GetEmbeddedAtoms()) {
      if (eri1.first == eri2.first) {
        res.push_back(eri1.first);
      }
    }
  }
  return res;
}

void EmbeddedFrag::mergeNoCommon(EmbeddedFrag &embObj, unsigned int toAid,
                                 unsigned int nbrAid) {
  // merge embObj to this fragment when there are no common atoms between the
  // two fragments
  PRECONDITION(dp_mol, "");
  // check that both this fragment and the one we are merging with belong to the
  // same molecule
  PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
  RDKit::INT_VECT commAtms;
  this->addNonRingAtom(nbrAid, toAid);
  embObj.addNonRingAtom(toAid, nbrAid);
  commAtms.push_back(toAid);
  commAtms.push_back(nbrAid);
  this->mergeWithCommon(embObj, commAtms);
}

void EmbeddedFrag::mergeWithCommon(EmbeddedFrag &embObj,
                                   RDKit::INT_VECT &commAtms) {
  PRECONDITION(dp_mol, "");
  PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
  PRECONDITION(commAtms.size() >= 1, "");

  // we already have one or more common atoms between this fragment One atom in
  // common can happen (look at issue 173)
  // - for cases where a cis/trans double bond is being merged with a ring
  //   system that shares one of atoms on the double bond.
  // - or if 'this' fragment was created by user specified coordinates - where
  //   only part of a fused ring system or cis/trans system was specified

  // if we have one atom in common, we have to deal with it carefully -
  unsigned int ctCase =
      0;  // book-keeper - if we have to merge a ring with a cis/trans dbl bond
  // what kind is it
  //    0 - if we are doing a cis/trans merge,
  //    1 - cis/trans and embObj is the dblBond,
  //    2 - cis/trans merge and 'this' is the dblBond

  if (commAtms.size() == 1) {
    // couple of possibilities here
    // 1. we are merging a ring system with a cis/trans dbl bond
    // 2. We are merging with a fused ring system out of which one of the atoms
    //    has already been embedded because the user specified its coordinates
    // First deal with the cis/trans case
    auto commAid = commAtms.front();
    int otherAtom = -1;
    if (d_eatoms[commAid].CisTransNbr >= 0) {
      ctCase = 2;
      // this fragment is the cis/trans dbl bond
      otherAtom =
          d_eatoms[commAid].nbr1;  // this is the other atom on the double bnd
      // now add this atom to the other fragment
      embObj.addNonRingAtom(otherAtom, commAid);  //, mol);
    } else if (embObj.d_eatoms[commAid].CisTransNbr >= 0) {
      ctCase = 1;
      // otherwise embObj is the cis/trans dbl bond
      otherAtom = embObj.d_eatoms[commAid].nbr1;
      this->addNonRingAtom(otherAtom, commAid);  //, mol);
    } else {
      otherAtom = d_eatoms[commAid].nbr1;
      if (otherAtom >= 0) {
        embObj.addNonRingAtom(otherAtom, commAid);  //, mol);
      }
    }
    if (otherAtom >= 0) {
      commAtms.push_back(otherAtom);
    }
  }

  RDGeom::Transform2D rtrans;
  if (commAtms.size() == 1) {
    // if we have only one atom in common we will use a one atom transform
    rtrans.assign(this->computeOneAtomTrans(commAtms.front(), embObj));
  } else {
    // if we have more than one we will use a two point transform
    auto cid1 = commAtms[0];
    auto cid2 = commAtms[1];
    const auto &ref1 = d_eatoms.at(cid1).loc;
    const auto &ref2 = d_eatoms.at(cid2).loc;
    const auto &oth1 = embObj.GetEmbeddedAtom(cid1).loc;
    const auto &oth2 = embObj.GetEmbeddedAtom(cid2).loc;
    // now compute the transform
    rtrans.SetTransform(ref1, ref2, oth1, oth2);
  }

  // transform the second fragment
  embObj.Transform(rtrans);

  // check to see if this transform screws up any cis/trans specifications
  if (commAtms.size() >= 2) {
    if (ctCase > 0) {
      // we have a cis/trans case we may have violated the specification
      // check and correct it with a reflection
      reflectIfNecessaryCisTrans(embObj, ctCase, commAtms[0], commAtms[1]);
    } else if (commAtms.size() == 2) {
      // we have just two atoms in common but we may a simply overcrowed one
      // side check for crowding and reflect
      reflectIfNecessaryDensity(embObj, commAtms[0], commAtms[1]);
    } else {
      // finally if we have more than two atoms in common - we will use the
      // third atom to figure out if we need a reflection12
      reflectIfNecessaryThirdPt(embObj, commAtms[0], commAtms[1], commAtms[2]);
    }
  }

  // finally merge the fragment by copying the non common atoms
  const auto &oatoms = embObj.GetEmbeddedAtoms();
  // copy the eatoms in embObj to this fragment
  for (const auto &ori : oatoms) {
    auto aid = ori.first;
    if (std::find(commAtms.begin(), commAtms.end(), aid) == commAtms.end()) {
      d_eatoms[aid] = ori.second;
      // also if any of these atoms have unattached neighbors add them to the
      // queue
      if (!ori.second.neighs.empty()) {
        if (std::find(d_attachPts.begin(), d_attachPts.end(), aid) ==
            d_attachPts.end()) {
          d_attachPts.push_back(aid);
        }
      }
    } else {
      if (ori.second.CisTransNbr >= 0) {
        d_eatoms[aid].CisTransNbr = ori.second.CisTransNbr;
        d_eatoms[aid].normal = ori.second.normal;
        d_eatoms[aid].ccw = ori.second.ccw;
      }
      if (ori.second.angle > 0.0) {
        d_eatoms[aid].angle = ori.second.angle;

        d_eatoms[aid].nbr1 = ori.second.nbr1;
        d_eatoms[aid].nbr2 = ori.second.nbr2;
      }
    }
  }

  // remember to update the not yet done neighbor of nbrAid
  for (auto cai : commAtms) {
    this->updateNewNeighs(cai);
  }
}

void EmbeddedFrag::mergeFragsWithComm(std::list<EmbeddedFrag> &efrags) {
  PRECONDITION(dp_mol, "");
  // first merge any fragments what share atoms in common
  auto nfri = efrags.end();
  while (1) {
    RDKit::INT_VECT commAtms;
    for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
      if (!efri->isDone()) {
        commAtms = this->findCommonAtoms(*efri);
        if (commAtms.size() > 0) {
          nfri = efri;
          break;
        }
      }
    }
    if (commAtms.empty()) {
      break;
    }

    CHECK_INVARIANT(nfri != efrags.end(), "iterator not initialized");
    this->mergeWithCommon((*nfri), commAtms);  //, mol);
    for (auto cai : commAtms) {
      if (d_eatoms.at(cai).neighs.empty() &&
          (std::find(d_attachPts.begin(), d_attachPts.end(), cai) !=
           d_attachPts.end())) {
        d_attachPts.erase(
            std::remove(d_attachPts.begin(), d_attachPts.end(), cai));
      }
    }
    efrags.erase(nfri);
  }
}

void EmbeddedFrag::expandEfrag(RDKit::INT_LIST &nratms,
                               std::list<EmbeddedFrag> &efrags) {
  PRECONDITION(dp_mol, "");

  // first merge any fragments that share atoms in common

  this->mergeFragsWithComm(efrags);  //, dp_mol);

  while (d_attachPts.size() > 0) {
    auto aid = d_attachPts.front();
    auto nbrs = d_eatoms[aid].neighs;
    CHECK_INVARIANT(!nbrs.empty(), "");
    for (auto nbri : nbrs) {
      auto nratmi = std::find(nratms.begin(), nratms.end(), nbri);
      if (nratmi != nratms.end()) {
        // the neighbor we have to add is a non ring atoms
        this->addNonRingAtom(nbri, aid);  //, mol);
        // remove this atom we just added from the nnratms list
        nratms.erase(nratmi);
      } else {
        // the neighbor atom must be part of a different embedded fragment -
        // merge that fragment with this one
        auto nfri = efrags.end();
        for (auto efri = efrags.begin(); efri != efrags.end(); ++efri) {
          // don't search fragments that are done
          if (!efri->isDone()) {
            const auto &eatoms = efri->GetEmbeddedAtoms();
            if (eatoms.find(nbri) != eatoms.end()) {
              nfri = efri;
              break;
            }
          }
        }
        if (nfri != efrags.end()) {
          this->mergeNoCommon((*nfri), aid, nbri);  //, mol);
          if (d_eatoms.at(nbri).neighs.empty() &&
              (std::find(d_attachPts.begin(), d_attachPts.end(), nbri) !=
               d_attachPts.end())) {
            d_attachPts.erase(
                std::remove(d_attachPts.begin(), d_attachPts.end(), nbri));
          }
          // remove this fragment from the list of embedded fragments
          efrags.erase(nfri);
        }
      }
    }

    // ok we are done with this atom forever
    d_attachPts.pop_front();
    d_eatoms[aid].neighs.clear();
    // now that we added new atoms to the this fragments - check if there are
    // new fragment we have common atoms with and merge with them
    this->mergeFragsWithComm(efrags);  //, mol);
  }
}

void EmbeddedFrag::Transform(const RDGeom::Transform2D &trans) {
  for (auto &eri : d_eatoms) {
    eri.second.Transform(trans);
  }
}

void EmbeddedFrag::computeBox() {
  d_px = -1.0e8;
  d_nx = 1.0e8;
  d_py = -1.0e8;
  d_ny = 1.0e8;

  for (const auto &eri : d_eatoms) {
    const auto &loc = eri.second.loc;
    d_px = std::max(d_px, loc.x);
    d_nx = std::min(d_nx, loc.x);
    d_py = std::max(d_py, loc.y);
    d_ny = std::min(d_ny, loc.y);
  }
  d_nx *= -1.0;
  d_ny *= -1.0;
}

void EmbeddedFrag::canonicalizeOrientation() {
  // fix for issue 198
  // no need to canonicalize if we are dealing with a single atm
  if (d_eatoms.size() <= 1) {
    return;
  }

  RDGeom::Point2D cent(0.0, 0.0);
  for (const auto &elem : d_eatoms) {
    cent += elem.second.loc;
  }
  cent *= (1.0 / d_eatoms.size());

  double xx = 0.0;
  double xy = 0.0;
  double yy = 0.0;

  // shift the center of the fragment to the origin and compute the covariance
  // matrix
  for (auto &elem : d_eatoms) {
    elem.second.loc -= cent;
    xx += (elem.second.loc.x) * (elem.second.loc.x);
    xy += (elem.second.loc.x) * (elem.second.loc.y);
    yy += (elem.second.loc.y) * (elem.second.loc.y);
  }

  RDGeom::Point2D eig1, eig2;
  // the eigen vectors are given by
  //   (2*xy, (yy - xx) + d) and (2*xy, (yy - xx) - d)
  // where d = sqrt((xx - yy)^2 + 4*xy^2)
  auto d = (xx - yy) * (xx - yy) + 4 * xy * xy;
  d = sqrt(d);
  eig1.x = 2 * xy;
  eig1.y = (yy - xx) + d;
  if (eig1.length() <= 1e-4) {
    return;
  }
  auto eVal1 = (xx + yy + d) / 2;
  eig1.normalize();

  eig2.x = 2 * xy;
  eig2.y = (yy - xx) - d;
  auto eVal2 = (xx + yy - d) / 2;

  if (eig2.length() > 1e-4) {
    eig2.normalize();

    // make sure eig1 corresponds to the larger eigenvalue:
    if (eVal2 > eVal1) {
      std::swap(eig1, eig2);
    }
  }
  // now rotate eig1 onto the X axis:
  RDGeom::Transform2D trans;
  trans.setVal(0, 0, eig1.x);
  trans.setVal(1, 0, -eig1.y);
  trans.setVal(0, 1, eig1.y);
  trans.setVal(1, 1, eig1.x);
  this->Transform(trans);
}

void _recurseAtomOneSide(unsigned int endAid, unsigned int begAid,
                         const RDKit::ROMol *mol, RDKit::INT_VECT &flipAids) {
  PRECONDITION(mol, "");
  flipAids.push_back(endAid);
  for (auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(endAid))) {
    if (nbr->getIdx() != begAid &&
        (std::find(flipAids.begin(), flipAids.end(),
                   static_cast<int>(nbr->getIdx())) == flipAids.end())) {
      _recurseAtomOneSide(nbr->getIdx(), begAid, mol, flipAids);
    }
  }
  return;
}

double _crossVal(const RDGeom::Point2D &v1, const RDGeom::Point2D &v2) {
  return v1.x * v2.y - v2.x * v1.y;
}

int _pairDIICompAscending(const PAIR_D_I_I &arg1, const PAIR_D_I_I &arg2) {
  return (arg1.first < arg2.first);
}

PAIR_I_I _findClosestPair(unsigned int beg1, unsigned int end1,
                          unsigned int beg2, unsigned int end2,
                          const RDKit::ROMol &mol, const double *dmat) {
  auto na = mol.getNumAtoms();
  auto d1 = dmat[beg1 * na + beg2];
  auto d2 = dmat[beg1 * na + end2];
  auto d3 = dmat[end1 * na + beg2];
  auto d4 = dmat[end1 * na + end2];
  auto minPr =
      std::min(PAIR_D_I_I(d1, PAIR_I_I(beg1, beg2)),
               PAIR_D_I_I(d2, PAIR_I_I(beg1, end2)), _pairDIICompAscending);
  minPr = std::min(minPr, PAIR_D_I_I(d3, PAIR_I_I(end1, beg2)),
                   _pairDIICompAscending);
  minPr = std::min(minPr, PAIR_D_I_I(d4, PAIR_I_I(end1, end2)),
                   _pairDIICompAscending);
  return minPr.second;
}

void EmbeddedFrag::computeDistMat(DOUBLE_SMART_PTR &dmat) {
  auto dmatPtr = dmat.get();
  for (auto efi = d_eatoms.begin(); efi != d_eatoms.end(); ++efi) {
    auto pti = efi->second.loc;
    auto ai = efi->first;
    for (auto efj = d_eatoms.begin(); efj != efi; ++efj) {
      auto ptj = efj->second.loc;
      auto aj = efj->first;
      ptj -= pti;
      if (ai < aj) {
        std::swap(ai, aj);
      }
      dmatPtr[(ai * (ai - 1) / 2) + aj] = ptj.length();
    }
  }
}

double EmbeddedFrag::mimicDistMatAndDensityCostFunc(
    const DOUBLE_SMART_PTR *dmat, double mimicDmatWt) {
  const double *ddata;
  if (dmat) {
    ddata = dmat->get();
  } else {
    ddata = nullptr;
  }
  auto na = dp_mol->getNumAtoms();
  if (na < 2) {
    return 0;
  }
  auto dsize = na * (na - 1) / 2;
  auto *ddata2D = new double[dsize];
  DOUBLE_SMART_PTR dmat2D(ddata2D);
  this->computeDistMat(dmat2D);
  double res1 = 0.0;
  double res2 = 0.0;
  for (auto i = 0u; i < dsize; ++i) {
    auto d = ddata2D[i];
    auto d2 = d * d;
    if (d2 > 1.e-3) {
      res1 += 1.0 / d2;
    } else {
      res1 += 1000.0;
    }
    if (ddata && (ddata[i] >= 0.0)) {
      auto dd = d - ddata[i];
      res2 += dd * dd;
    }
  }

  auto wt = mimicDmatWt;
  if (wt > 1.0) {
    wt = 1.0;
  } else if (wt < 0.0) {
    wt = 0.0;
  }

  return ((1.0 - wt) * res1) + (wt * res2);
}

// Permute the bonds at a degree 4 node
//
//      A                    B
//      |                    |
//   B--C--D     to       A--C--D
//      |                    |
//      E                    E
//
// Note that everything attached to B and A are also effected. This is what
// happens here
// 1. Find the line "l" bisecting the angle BCA
// 2. Find the atoms in the fragment generated by breaking the bond between C
//    and A that includes A. Lets call is Fa
// 3. Similarly find the fragment Fb that includes B by breaking the bond CB
// 4. Reflect Fb and Fa through "l"
void EmbeddedFrag::permuteBonds(unsigned int aid, unsigned int aid1,
                                unsigned int aid2) {
  PRECONDITION(dp_mol, "");
  auto rl1 = d_eatoms.at(aid).loc;
  auto rl2 = d_eatoms.at(aid1).loc + d_eatoms.at(aid2).loc;
  rl2 *= 0.5;

  RDKit::INT_VECT fragA, fragB;

  // now find the fragment that contains aid1 but not aid
  _recurseAtomOneSide(aid1, aid, dp_mol, fragA);

  // now find the fragment that contains aid2 but not aid
  _recurseAtomOneSide(aid2, aid, dp_mol, fragB);

  // now just loop through these atoms and reflect them
  for (auto fi : fragA) {
    d_eatoms[fi].Reflect(rl1, rl2);
  }

  for (auto fi : fragB) {
    d_eatoms[fi].Reflect(rl1, rl2);
  }
}

void EmbeddedFrag::randomSampleFlipsAndPermutations(
    unsigned int nBondsPerSample, unsigned int nSamples, int seed,
    const DOUBLE_SMART_PTR *dmat, double mimicDmatWt, bool permuteDeg4Nodes) {
  PRECONDITION(dp_mol, "");

  const auto &rotBonds = getAllRotatableBonds(*dp_mol);
  auto nb = rotBonds.size();  // number of rotatable bonds that can be flipped

  // if we also want to permute deg 4 nodes, find out how many of these are
  // around and can be permuted
  unsigned int nd4 = 0;
  RDKit::INT_VECT deg4nodes;
  RDKit::VECT_INT_VECT deg4NbrBids, deg4NbrAids;

  if (permuteDeg4Nodes) {
    for (const auto atom : dp_mol->atoms()) {
      auto caid = atom->getIdx();
      if ((getDepictDegree(atom) == 4) &&
          (!(dp_mol->getRingInfo()->numAtomRings(caid)))) {
        RDKit::INT_VECT aids, bids;
        getNbrAtomAndBondIds(caid, dp_mol, aids, bids);
        // make sure all the atoms in aids are in this embeddedfrag and can be
        // perturbed
        bool allin = true;
        for (auto aid : aids) {
          auto nbrIter = d_eatoms.find(aid);
          if (nbrIter == d_eatoms.end() || nbrIter->second.df_fixed) {
            allin = false;
            break;
          }
        }
        if (allin) {
          deg4nodes.push_back(caid);
          deg4NbrBids.push_back(bids);
          deg4NbrAids.push_back(aids);
        }
      }
    }
    nd4 = deg4nodes.size();
  }

  unsigned int nt = nb + nd4;

  unsigned int nPerSample = std::min(nt, nBondsPerSample);

  auto &generator = RDKit::getRandomGenerator();
  if (seed > 0) {
    generator.seed(seed);
  }
  RDKit::uniform_int dist(0, nt - 1);
  RDKit::int_source_type intRandomSrc(generator, dist);

  RDGeom::INT_POINT2D_MAP bestCrdMap;
  auto bestDens = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt);
  for (const auto &efi : d_eatoms) {
    bestCrdMap[efi.first] = efi.second.loc;
  }
  for (auto si = 0u; si < nSamples; ++si) {
    // randomly pick nPerSample bonds and flip them
    for (auto fi = 0u; fi < nPerSample; ++fi) {
      unsigned int ri = intRandomSrc();
      // if ri is less than the number of rotatable bonds (nb), we will flip a
      // rot bond
      if (ri < nb) {
        this->flipAboutBond(rotBonds.at(ri));
      } else {  // ri is >= nb we permute the bonds at a deg 4 node
        unsigned int d4i =
            ri - nb;  // so we will permute at the 'di'th degree 4 node
        auto ai = deg4nodes.at(d4i);
        // collect the locations for the neighbors
        VECT_C_POINT nbrLocs;
        for (auto aci : deg4NbrAids[d4i]) {
          nbrLocs.push_back(&(d_eatoms.at(aci).loc));
        }
        auto bndPairs = findBondsPairsToPermuteDeg4(
            d_eatoms.at(ai).loc, deg4NbrBids.at(d4i), nbrLocs);

        auto rval = RDKit::getRandomVal();
        unsigned int fbi = 0;
        if (rval > 0.5) {
          fbi = 1;
        }
        auto aid1 =
            dp_mol->getBondWithIdx(bndPairs.at(fbi).first)->getOtherAtomIdx(ai);
        auto aid2 = dp_mol->getBondWithIdx(bndPairs.at(fbi).second)
                        ->getOtherAtomIdx(ai);
        this->permuteBonds(ai, aid1, aid2);
      }
    }

    // compute the density of the structure and check if it improved
    auto density = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt);
    // if (density < bestDens) {
    if (bestDens - density > 1e-4) {
      bestDens = density;
      for (const auto &efi : d_eatoms) {
        bestCrdMap[efi.first] = efi.second.loc;
      }
    }
  }
  // now copy the best coordinates to the fragment
  for (auto &efi : d_eatoms) {
    efi.second.loc = bestCrdMap.at(efi.first);
  }
}

std::vector<PAIR_I_I> EmbeddedFrag::findCollisions(const double *dmat,
                                                   bool includeBonds) {
  // find a pair of atoms that are too close to each other
  std::vector<PAIR_I_I> res;
  for (auto &d_eatom : d_eatoms) {
    d_eatom.second.d_density = 0.0;
  }

  auto colThres2 = COLLISION_THRES * COLLISION_THRES;
  // if we a re dealing with non carbon atoms we will increase the collision
  // threshold. This is because only hetero atoms are typically drawn in a
  // depiction.
  double atomTypeFactor1, atomTypeFactor2;
  for (auto efi = d_eatoms.begin(); efi != d_eatoms.end(); ++efi) {
    atomTypeFactor1 = 1.0;
    if (dp_mol->getAtomWithIdx(efi->first)->getAtomicNum() != 6) {
      atomTypeFactor1 = HETEROATOM_COLL_SCALE;
    }
    for (auto efj = d_eatoms.begin(); efj != efi; ++efj) {
      atomTypeFactor2 = 1.0;
      if (dp_mol->getAtomWithIdx(efj->first)->getAtomicNum() != 6) {
        atomTypeFactor2 = HETEROATOM_COLL_SCALE;
      }
      auto ptj = efj->second.loc;
      ptj -= efi->second.loc;
      auto d2 = ptj.lengthSq();
      if (d2 > 1.0e-3) {
        efi->second.d_density += (1 / d2);
        efj->second.d_density += (1 / d2);
      } else {
        efi->second.d_density += 1000.0;
        efj->second.d_density += 1000.0;
      }
      d2 /= (atomTypeFactor1 * atomTypeFactor2);
      if (d2 < colThres2) {
        PAIR_I_I cAids(efi->first, efj->first);
        res.push_back(cAids);
      }
    }
  }
  if (includeBonds) {
    // now find bond collisions
    double BOND_THRES2 = BOND_THRES * BOND_THRES;
    for (const auto b1 : dp_mol->bonds()) {
      auto bid1 = b1->getIdx();
      auto beg1 = b1->getBeginAtomIdx();
      auto end1 = b1->getEndAtomIdx();
      if ((d_eatoms.find(beg1) != d_eatoms.end()) &&
          (d_eatoms.find(end1) != d_eatoms.end())) {
        auto v1 = d_eatoms[end1].loc - d_eatoms[beg1].loc;
        auto avg1 = d_eatoms[end1].loc + d_eatoms[beg1].loc;
        avg1 *= 0.5;
        for (const auto b2 : dp_mol->bonds()) {
          if (b2->getIdx() <= bid1) {
            continue;
          }

          auto beg2 = b2->getBeginAtomIdx();
          auto end2 = b2->getEndAtomIdx();
          if ((d_eatoms.find(beg2) != d_eatoms.end()) &&
              (d_eatoms.find(end2) != d_eatoms.end())) {
            auto avg2 = d_eatoms[end2].loc + d_eatoms[beg2].loc;
            avg2 *= 0.5;
            avg2 -= avg1;
            if (avg2.lengthSq() < 0.5 && avg2.lengthSq() < BOND_THRES2) {
              auto v2 = d_eatoms[beg2].loc - d_eatoms[beg1].loc;
              auto v3 = d_eatoms[end2].loc - d_eatoms[beg1].loc;
              auto valProd = _crossVal(v1, v2) * _crossVal(v1, v3);
              if (valProd < -1e-6) {
                // we have a collision, find the closest two atoms
                auto cAids =
                    _findClosestPair(beg1, end1, beg2, end2, *dp_mol, dmat);
                res.push_back(cAids);
              }
            }
          }
        }
      }
    }
  }
  return res;
}

double EmbeddedFrag::totalDensity() {
  return std::accumulate(
      d_eatoms.begin(), d_eatoms.end(), 0.0,
      [](double accum, auto &dea) { return dea.second.d_density + accum; });
}

void _recurseDegTwoRingAtoms(unsigned int aid, const RDKit::ROMol *mol,
                             RDKit::INT_VECT &rPath,
                             RDKit::INT_INT_VECT_MAP &nbrMap) {
  PRECONDITION(mol, "");
  // find all atoms along a path that have two ring atoms on them
  // aid is where will start looking and then we will recurse
  RDKit::INT_VECT nbrs;
  for (const auto bnd : mol->atomBonds(mol->getAtomWithIdx(aid))) {
    if (mol->getRingInfo()->numBondRings(bnd->getIdx())) {
      nbrs.push_back(bnd->getOtherAtomIdx(aid));
    }
  }
  if (nbrs.size() != 2) {
    return;
  } else {
    rPath.push_back(aid);
    nbrMap[aid] = nbrs;
    for (auto nbr : nbrs) {
      if (std::find(rPath.begin(), rPath.end(), nbr) == rPath.end()) {
        _recurseDegTwoRingAtoms(nbr, mol, rPath, nbrMap);
      }
    }
  }
}

unsigned int _anyNonRingBonds(unsigned int aid, RDKit::INT_LIST path,
                              const RDKit::ROMol *mol) {
  PRECONDITION(mol, "");
  // check if there are any non-ring bonds on the path starting at aid
  auto prev = aid;
  auto nOpen = 0u;
  for (auto pi : path) {
    const auto bond = mol->getBondBetweenAtoms(prev, pi);
    CHECK_INVARIANT(bond, "no bond found");
    if (!mol->getRingInfo()->numBondRings(bond->getIdx())) {
      ++nOpen;
    }
    prev = pi;
  }
  return nOpen;
}

void EmbeddedFrag::flipAboutBond(unsigned int bondId, bool flipEnd) {
  PRECONDITION(dp_mol, "");
  PRECONDITION(bondId < dp_mol->getNumBonds(), "");
  // reflect all the atoms on one side of a bond using the bond as the mirror
  const auto bond = dp_mol->getBondWithIdx(bondId);

  // we should not be flip things around a ring bond
  CHECK_INVARIANT(!(dp_mol->getRingInfo()->numBondRings(bondId)), "");

  auto begAid = bond->getBeginAtomIdx();
  auto endAid = bond->getEndAtomIdx();

  if (!flipEnd) {
    std::swap(begAid, endAid);
  }

  const auto &begLoc = d_eatoms.at(begAid).loc;
  const auto &endLoc = d_eatoms.at(endAid).loc;

  // arbitrary choice here - find all atoms on one side of the bond
  // endAtom side - we will do this recursively
  RDKit::INT_VECT endSideAids;
  _recurseAtomOneSide(endAid, begAid, dp_mol, endSideAids);

  // look for fixed atoms in the fragment:
  unsigned int nAtomsFixed = 0;
  for (auto &d_eatom : d_eatoms) {
    if (d_eatom.second.df_fixed) {
      ++nAtomsFixed;
    }
  }
  // if there are fixed atoms, look at the atoms on the "end side"
  unsigned int nEndAtomsFixed = 0;
  if (nAtomsFixed) {
    for (auto endAtomId : endSideAids) {
      if (d_eatoms[endAtomId].df_fixed) {
        ++nEndAtomsFixed;
      }
    }
  }
  // now we have the molecule split into two groups of atoms
  // atom on the side of endAid and the rest.
  // we will flip the side that is smaller, assuming that there
  // are no fixed atoms there
  bool endSideFlip = true;
  if (nEndAtomsFixed) {
    endSideFlip = false;
    // there are fixed atoms on both sides, just return
    return;
  } else {
    auto nats = d_eatoms.size();
    auto nEndSide = endSideAids.size();
    if ((nats - nEndSide) < nEndSide) {
      endSideFlip = false;
    }
  }
  for (auto &d_eatom : d_eatoms) {
    const auto fii = std::find(endSideAids.begin(), endSideAids.end(),
                               static_cast<int>(d_eatom.first));
    if (endSideFlip ^ (fii == endSideAids.end())) {
      d_eatom.second.Reflect(begLoc, endLoc);
    }
  }
}

unsigned int _findDeg1Neighbor(const RDKit::ROMol *mol, unsigned int aid) {
  PRECONDITION(mol, "");
  auto deg = getDepictDegree(mol->getAtomWithIdx(aid));
  CHECK_INVARIANT(deg == 1, "");
  return *mol->getAtomNeighbors(mol->getAtomWithIdx(aid)).first;
}

unsigned int _findClosestNeighbor(const RDKit::ROMol *mol, const double *dmat,
                                  unsigned int aid1, unsigned int aid2) {
  PRECONDITION(mol, "");
  unsigned int res = 0;
  double mdist = 1.e8;
  auto naid = aid1 * (mol->getNumAtoms());
  for (const auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(aid2))) {
    auto d = dmat[naid + nbr->getIdx()];
    if (d < mdist) {
      mdist = d;
      res = nbr->getIdx();
    }
  }
  return res;
}

void EmbeddedFrag::openAngles(const double *dmat, unsigned int aid1,
                              unsigned int aid2) {
  // Assuming that either aid1, and/or aid2 are degree 1 atoms, we will open up
  // the angles
  //
  //     1 2
  //    /   \                                                   this space
  //   /     \                                       intentionally left blank
  //  a-------b
  //
  // If 1 and 2 are too close to each other we open up angle(1ab) if 1 is a
  // degree 1 node and
  // angle(2ba) if 2 is a degree 1 node. Say 1 is a degree 1 node but 2 is not.
  // Then from the neighbors of 2 we need to choose which one should be b. Also
  // keep in mind
  // that a need not be a neighbor of b. In this case we will pick b to be the
  // closest neighbor of a

  PRECONDITION(dp_mol, "");
  PRECONDITION(dmat, "");
  auto deg1 = getDepictDegree(dp_mol->getAtomWithIdx(aid1));
  auto deg2 = getDepictDegree(dp_mol->getAtomWithIdx(aid2));
  auto fixed1 = d_eatoms.at(aid1).df_fixed;
  auto fixed2 = d_eatoms.at(aid2).df_fixed;
  if ((deg1 > 1 || fixed1) && (deg2 > 1 || fixed2)) {
    return;
  }
  unsigned int aidA;
  unsigned int aidB;
  int type = 0;
  if ((deg1 == 1 && !fixed1) && (deg2 == 1 && !fixed2)) {
    aidA = _findDeg1Neighbor(dp_mol, aid1);
    aidB = _findDeg1Neighbor(dp_mol, aid2);
    type = 1;
  } else if ((deg1 == 1 && !fixed1) && (deg2 > 1 || fixed2)) {
    aidA = _findDeg1Neighbor(dp_mol, aid1);
    aidB = _findClosestNeighbor(dp_mol, dmat, aidA, aid2);
    type = 2;
  } else {
    aidB = _findDeg1Neighbor(dp_mol, aid2);
    aidA = _findClosestNeighbor(dp_mol, dmat, aidB, aid1);
    type = 3;
  }

  auto v2 = d_eatoms.at(aid1).loc - d_eatoms.at(aidA).loc;
  auto v1 = d_eatoms.at(aidB).loc - d_eatoms.at(aidA).loc;
  auto cross = (v1.x) * (v2.y) - (v1.y) * (v2.x);
  double angle;
  RDGeom::Transform2D trans1, trans2;
  switch (type) {
    case 1:
      angle = ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans1.SetTransform(d_eatoms[aidA].loc, angle);
      trans2.SetTransform(d_eatoms[aidB].loc, -1.0 * angle);
      trans1.TransformPoint(d_eatoms[aid1].loc);
      trans2.TransformPoint(d_eatoms[aid2].loc);
      break;
    case 2:
      angle = 2.0 * ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans1.SetTransform(d_eatoms[aidA].loc, angle);
      trans1.TransformPoint(d_eatoms[aid1].loc);
      break;
    case 3:
      angle = -2.0 * ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans2.SetTransform(d_eatoms[aidB].loc, angle);
      trans2.TransformPoint(d_eatoms[aid2].loc);
      break;
    default:
      break;
  }
}

void EmbeddedFrag::removeCollisionsBondFlip() {
  // try to remove collisions in a structure by flipping rotatable bonds along
  // the shortest path between the colliding atoms. we will limit the number of
  // times we are going to do this since we may fall into spiral where removing
  // a collision may create a new one
  auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
  auto colls = this->findCollisions(dmat);
  std::map<int, unsigned int> doneBonds;
  unsigned int iter = 0;
  while (iter < MAX_COLL_ITERS && colls.size()) {
    auto ncols = colls.size();
    if (ncols > 0) {
      // we have a collision
      auto cAids = colls[0];
      auto rotBonds = getRotatableBonds(*dp_mol, cAids.first, cAids.second);
      auto prevDensity = this->totalDensity();
      for (auto ri : rotBonds) {
        auto doneBondsRiIt = doneBonds.find(ri);
        if ((doneBondsRiIt == doneBonds.end()) ||
            (doneBondsRiIt->second < NUM_BONDS_FLIPS)) {
          if (doneBondsRiIt == doneBonds.end()) {
            doneBonds[ri] = 1;
          } else {
            doneBondsRiIt->second += 1;
          }

          flipAboutBond(ri);
          colls = this->findCollisions(dmat);
          auto newDensity = this->totalDensity();
          if (colls.size() < ncols) {
            doneBonds[ri] = NUM_BONDS_FLIPS;  // lock this rotatable bond
            break;
          } else if (colls.size() == ncols && newDensity < prevDensity) {
            break;
          } else {
            // we made the wrong move earlier - reject the flip move it back
            flipAboutBond(ri);
            colls = this->findCollisions(dmat);
            // and try the other end:
            flipAboutBond(ri, false);
            colls = this->findCollisions(dmat);
            newDensity = this->totalDensity();
            if (colls.size() < ncols) {
              doneBonds[ri] = NUM_BONDS_FLIPS;  // lock this rotatable bond
              break;
            } else if (colls.size() == ncols && newDensity < prevDensity) {
              break;
            } else {
              flipAboutBond(ri, false);
              colls = this->findCollisions(dmat);
            }
          }
        }
      }
    }
    ++iter;
  }
}

void EmbeddedFrag::removeCollisionsOpenAngles() {
  auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
  // try opening up angles
  for (const auto &cpi : this->findCollisions(dmat, 0)) {
    // find out which of the two offending atoms we want to move
    // we will use the one with the smallest degree
    this->openAngles(dmat, cpi.first, cpi.second);
  }
}

void EmbeddedFrag::removeCollisionsShortenBonds() {
  auto dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
  // if there are still some collision points left - flipping rotatable bonds
  // and opening angles is not doing it - we will try two last things
  //  - if all the bonds between the colliding atoms are rings bonds,
  //    we most likely have a collision within a bridged system (Issue 199).
  //    In this case we will try to find a path of colliding atoms (in one
  //    of the rings) and shorten all the bond in the path
  //  - on the other hand if we have non-ring bonds as well in the path
  //    between the colliding atoms we will simply shorten each one of
  //    them by a little bit.
  auto colls = this->findCollisions(dmat, 0);
  auto ncols = colls.size();
  auto iter = 0u;
  while (ncols && iter < MAX_COLL_ITERS) {
    const auto cAids = colls.front();
    // find out which of the two offending atoms we want to move
    // we will use the one with the smallest degree
    auto aid1 = cAids.first;
    auto aid2 = cAids.second;
    auto fixed1 = d_eatoms.at(aid1).df_fixed;
    auto fixed2 = d_eatoms.at(aid2).df_fixed;
    if (fixed1 && fixed2) {
      // both atoms are fixed, so there's nothing
      // we can do about this collision.
      colls.erase(colls.begin());
      ncols = colls.size();
      ++iter;
      continue;
    }
    auto deg1 = dp_mol->getAtomWithIdx(aid1)->getDegree();
    auto deg2 = dp_mol->getAtomWithIdx(aid2)->getDegree();
    if (fixed1 || (deg2 > deg1 && !fixed2)) {
      // reverse the order
      std::swap(deg1, deg2);
      std::swap(aid1, aid2);
      std::swap(fixed1, fixed2);
    }
    // now find the path between the two ends
    auto path = RDKit::MolOps::getShortestPath(*dp_mol, aid1, aid2);
    if (!path.size()) {
      // there's no path between the ends, so there's nothing
      // we can really do about this collision.
      colls.erase(colls.begin());
    } else {
      // aid1 is on the front of the path, pop it off:
      CHECK_INVARIANT(path.front() == aid1, "bad path head");
      path.pop_front();

      auto nOpen = _anyNonRingBonds(aid1, path, dp_mol);
      if (nOpen > 0) {
        if (deg1 == 1) {
          auto loc = d_eatoms.at(aid1).loc;
          auto aidA = _findDeg1Neighbor(dp_mol, aid1);
          loc -= d_eatoms[aidA].loc;
          loc *= .9;
          if (loc.length() > .75) {
            loc += d_eatoms[aidA].loc;
            d_eatoms[aid1].loc = loc;
          }
        }
        if (deg2 == 1 && !fixed2) {
          auto loc = d_eatoms.at(aid2).loc;
          auto aidA = _findDeg1Neighbor(dp_mol, aid2);
          loc -= d_eatoms[aidA].loc;
          loc *= .9;
          if (loc.length() > .75) {
            loc += d_eatoms[aidA].loc;
            d_eatoms[aid2].loc = loc;
          }
        }
      } else {
        // we probably have a bridged system
        // lets hope that aids has only two ring bond on it
        RDKit::INT_VECT rPath;
        RDKit::INT_INT_VECT_MAP nbrMap;
        _recurseDegTwoRingAtoms(aid1, dp_mol, rPath, nbrMap);
        if (rPath.size() == 0) {
          _recurseDegTwoRingAtoms(aid2, dp_mol, rPath, nbrMap);
        }
        // now we will take each of the atoms in rPath and
        // "move them in" a little bit this is what "move them
        //  in" means (what we need is hand drawn picture in the comments)
        // - let r1 and r2 be the ring neighbor of the current atom r0
        // - we will find the vector that bisects angle(r1, r0, r2)
        // - we will move r0 along this vector
        RDGeom::INT_POINT2D_MAP moveMap;
        for (auto rpi : rPath) {
          if (d_eatoms.at(rpi).df_fixed) {
            continue;
          }
          auto mv = d_eatoms[nbrMap[rpi][0]].loc;
          mv += d_eatoms[nbrMap.at(rpi)[1]].loc;
          mv *= 0.5;
          mv -= d_eatoms.at(rpi).loc;
          mv.normalize();
          mv *= COLLISION_THRES;
          moveMap[rpi] = mv;
        }
        for (auto rpi : rPath) {
          d_eatoms[rpi].loc += moveMap[rpi];
        }
      }
      colls = this->findCollisions(dmat, 0);
    }
    ncols = colls.size();
    ++iter;
  }
}
}  // namespace RDDepict
