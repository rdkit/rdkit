//
//  Copyright (C) 2015-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/SubstanceGroup.h>

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ****************************************************************************

namespace RDKit {
namespace MolDraw2D_detail {
// implementation from $RDBASE/rdkit/sping/pid.py
void arcPoints(const Point2D &cds1, const Point2D &cds2,
               std::vector<Point2D> &res, float startAng, float extent) {
  // Note: this implementation is simple and not particularly efficient.
  float xScale = (cds2.x - cds1.x) / 2.0;
  float yScale = (cds2.y - cds1.y) / 2.0;
  if (xScale < 0) {
    xScale *= -1;
  }
  if (yScale < 0) {
    yScale *= -1;
  }

  float x = std::min(cds1.x, cds2.x) + xScale;
  float y = std::min(cds1.y, cds2.y) + yScale;

  int steps = std::max(static_cast<int>(extent * 2), 5);
  float step = M_PI * extent / (180 * steps);
  float angle = M_PI * startAng / 180;
  for (int i = 0; i <= steps; ++i) {
    Point2D point(x + xScale * cos(angle), y - yScale * sin(angle));
    res.emplace_back(point);
    angle += step;
  }
}

void addStereoAnnotation(const ROMol &mol, bool includeRelativeCIP) {
  const auto &sgs = mol.getStereoGroups();
  std::vector<unsigned int> doneAts(mol.getNumAtoms(), 0);
  unsigned int grpid = 1;
  for (const auto &sg : sgs) {
    for (const auto atom : sg.getAtoms()) {
      if (doneAts[atom->getIdx()]) {
        BOOST_LOG(rdWarningLog) << "Warning: atom " << atom->getIdx()
                                << " is in more than one stereogroup. Only the "
                                   "label from the first group will be used."
                                << std::endl;
        continue;
      }
      std::string lab;
      std::string cip;
      if (includeRelativeCIP ||
          sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        atom->getPropIfPresent(common_properties::_CIPCode, cip);
      }
      switch (sg.getGroupType()) {
        case StereoGroupType::STEREO_ABSOLUTE:
          lab = "abs";
          break;
        case StereoGroupType::STEREO_OR:
          lab = (boost::format("or%d") % grpid).str();
          break;
        case StereoGroupType::STEREO_AND:
          lab = (boost::format("and%d") % grpid).str();
          break;
        default:
          break;
      }
      if (!lab.empty()) {
        doneAts[atom->getIdx()] = 1;
        if (!cip.empty()) {
          lab += " (" + cip + ")";
        }
        atom->setProp(common_properties::atomNote, lab);
      }
    }
    if (sg.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
      ++grpid;
    }
  }
  for (auto atom : mol.atoms()) {
    std::string cip;
    if (!doneAts[atom->getIdx()] &&
        atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
      std::string lab = "(" + cip + ")";
      atom->setProp(common_properties::atomNote, lab);
    }
  }
  for (auto bond : mol.bonds()) {
    std::string cip;
    if (!bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
      if (bond->getStereo() == Bond::STEREOE) {
        cip = "E";
      } else if (bond->getStereo() == Bond::STEREOZ) {
        cip = "Z";
      }
    }
    if (!cip.empty()) {
      std::string lab = "(" + cip + ")";
      bond->setProp(common_properties::bondNote, lab);
    }
  }
}

void drawBracketsForSGroup(MolDraw2D &drawer, const ROMol &mol,
                           const SubstanceGroup &sg, const Conformer &conf) {
  std::vector<SubstanceGroup::Bracket> brackets = sg.getBrackets();
  if (brackets.empty()) {
    return;
  }
  RDGeom::Point2D center{0, 0};
  double avgLength = 0;
  for (const auto &brk : brackets) {
    const RDGeom::Point2D p1{brk[0]};
    const RDGeom::Point2D p2{brk[1]};
    auto v = p2 - p1;
    center += p1 + v / 2;
    avgLength += v.length();
  }
  center /= brackets.size();
  avgLength /= brackets.size();

  const auto ocolor = drawer.colour();
  const auto olw = drawer.lineWidth();
  const auto ffs = drawer.fontSize();
  drawer.setColour({0.4, 0.4, 0.4});
  drawer.setLineWidth(2);
  unsigned int whichBrk = 0;
  for (const auto &brk : brackets) {
    const RDGeom::Point2D p1{brk[0]};
    const RDGeom::Point2D p2{brk[1]};
    drawer.drawLine(p1, p2);
    auto v = p2 - p1;
    auto centerv = center - (p1 + v / 2);
    RDGeom::Point2D perp{v.y, -v.x};
    perp.normalize();
    if (perp.dotProduct(centerv) < 0) {
      perp *= -1;
    }
    perp *= avgLength / 10;
    drawer.drawLine(p1, p1 + perp);
    drawer.drawLine(p2, p2 + perp);
    ++whichBrk;
  }
  drawer.setColour(ocolor);
  drawer.setLineWidth(olw);
  drawer.setFontSize(ffs);
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
