//
//  Copyright (C) 2017-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <cstdlib>

#include "coordgen/sketcherMinimizer.h"
#include "coordgen/CoordgenFragmenter.h"

namespace RDKit {
namespace CoordGen {

struct CoordGenParams {
  const float sketcherCoarsePrecision = 0.01;
  const float sketcherStandardPrecision = SKETCHER_STANDARD_PRECISION;
  const float sketcherBestPrecision = SKETCHER_BEST_PRECISION;
  const float sketcherQuickPrecision = SKETCHER_QUICK_PRECISION;
  RDGeom::INT_POINT2D_MAP
      coordMap;  // coordinates for fixing particular atoms of a template
  const ROMol* templateMol = nullptr;  // a molecule to use as a template
  double coordgenScaling = 50.0;       // at the time this was written, coordgen
  // returned coordinates with a single bond
  // length of 50.
  std::string templateFileDir = "";
  float minimizerPrecision =
      sketcherCoarsePrecision;  // controls sketch precision
  bool dbg_useConstrained = true;   // debugging
  bool dbg_useFixed = false;        // debugging
  bool minimizeOnly = false;        // don't actually generate full coords
};

static CoordGenParams defaultParams;

//! Generates a 2D conformer for a molecule and replaces the existing conformers
/*!
  This call uses the CoordGen library from Schroedinger

  \param mol     the molecule we are working with
  \param params  a pointer to a parameter object

*/
template <typename T>
unsigned int addCoords(T& mol, const CoordGenParams* params = nullptr) {
  if (!params) params = &defaultParams;
  // FIX: the default value of this should be handled once in a threadsafe way
  std::string templateFileDir;
  if (params->templateFileDir != "") {
    templateFileDir = params->templateFileDir;
  } else {
    auto rdbase = std::getenv("RDBASE");
    if (rdbase != nullptr) {
      templateFileDir += rdbase;
      templateFileDir += "/Data/";
    }
  }

  if (params->minimizeOnly && !mol.getNumConformers()) {
    throw ValueErrorException(
        "minimizeOnly set but molecule has no conformers");
  }

  double scaleFactor = params->coordgenScaling;

  sketcherMinimizer minimizer(params->minimizerPrecision);
  auto min_mol = new sketcherMinimizerMolecule();

  // FIX: only do this check once.
  // std::cerr << "  TEMPLATES: " << templateFileDir << std::endl;
  if (templateFileDir != "") {
    minimizer.setTemplateFileDir(templateFileDir);
  }
  bool hasTemplateMatch = false;
  MatchVectType mv;
  if (!params->minimizeOnly && params->templateMol &&
      params->templateMol->getNumConformers() == 1) {
    if (SubstructMatch(mol, *(params->templateMol), mv)) {
      hasTemplateMatch = true;
    }
  }

  // if we're doing coordinate minimization it makes our life easier to
  // start by translating to the origin
  RDGeom::Point3D centroid{0.0, 0.0, 0.0};
  if (params->minimizeOnly) {
    auto conf = mol.getConformer();
    centroid = MolTransforms::computeCentroid(conf);
  }

  std::vector<sketcherMinimizerAtom*> ats(mol.getNumAtoms());
  for (auto atit = mol.beginAtoms(); atit != mol.endAtoms(); ++atit) {
    auto oatom = *atit;
    auto atom = min_mol->addNewAtom();
    atom->molecule = min_mol;  // seems like this should be in addNewAtom()
    atom->atomicNumber = oatom->getAtomicNum();
    atom->charge = oatom->getFormalCharge();
    if (params->minimizeOnly) {
      auto coords = mol.getConformer().getAtomPos(oatom->getIdx());
      coords -= centroid;
      atom->coordinates = sketcherMinimizerPointF(coords.x * scaleFactor,
                                                  coords.y * scaleFactor);
      atom->fixed = false;
      atom->constrained = false;
    } else if (hasTemplateMatch || params->coordMap.find(oatom->getIdx()) !=
                                       params->coordMap.end()) {
      atom->constrained = params->dbg_useConstrained;
      atom->fixed = params->dbg_useFixed;
      if (hasTemplateMatch) {
        for (auto& pr : mv) {
          if (pr.second == static_cast<int>(oatom->getIdx())) {
            const RDGeom::Point3D& coords =
                params->templateMol->getConformer().getAtomPos(pr.first);
            atom->templateCoordinates = sketcherMinimizerPointF(
                coords.x * scaleFactor, coords.y * scaleFactor);
            break;
          }
        }
      } else {
        const RDGeom::Point2D& coords =
            params->coordMap.find(oatom->getIdx())->second;
        atom->templateCoordinates = sketcherMinimizerPointF(
            coords.x * scaleFactor, coords.y * scaleFactor);
      }
    }
    ats[oatom->getIdx()] = atom;
  }

  // coordgen has its own ideas about what bond lengths should be
  // if we are doing minimizeOnly we want to do a rigid scaling after the
  // minimization in order to try and get as close to where we started as
  // possible
  double singleBondLenAccum = 0.0;
  double bondLenAccum = 0.0;
  unsigned nSingle = 0;

  std::vector<sketcherMinimizerBond*> bnds(mol.getNumBonds());
  for (auto bndit = mol.beginBonds(); bndit != mol.endBonds(); ++bndit) {
    auto obnd = *bndit;
    auto bnd = min_mol->addNewBond(ats[obnd->getBeginAtomIdx()],
                                   ats[obnd->getEndAtomIdx()]);
    // FIX: This is no doubt wrong
    switch (obnd->getBondType()) {
      case Bond::SINGLE:
        bnd->bondOrder = 1;
        if (params->minimizeOnly) {
          ++nSingle;
          singleBondLenAccum +=
              (mol.getConformer().getAtomPos(obnd->getBeginAtomIdx()) -
               mol.getConformer().getAtomPos(obnd->getEndAtomIdx()))
                  .length();
        }
        break;
      case Bond::DOUBLE:
        bnd->bondOrder = 2;
        break;
      case Bond::TRIPLE:
        bnd->bondOrder = 3;
        break;
      case Bond::AROMATIC:
        bnd->bondOrder = 1;
        break;
      default:
        BOOST_LOG(rdWarningLog) << "unrecognized bond type";
    }
    bnds[obnd->getIdx()] = bnd;
    if (params->minimizeOnly) {
      bondLenAccum += (mol.getConformer().getAtomPos(obnd->getBeginAtomIdx()) -
                       mol.getConformer().getAtomPos(obnd->getEndAtomIdx()))
                          .length();
    }
  }
  double avgBondLen = 1.0;
  if (params->minimizeOnly) {
    if (nSingle) {
      avgBondLen = singleBondLenAccum / nSingle;
    } else if (mol.getNumBonds()) {
      avgBondLen = bondLenAccum / mol.getNumBonds();
    }
  }

  // setup double bond stereo
  min_mol->assignBondsAndNeighbors(ats, bnds);
  for (auto bndit = mol.beginBonds(); bndit != mol.endBonds(); ++bndit) {
    auto obnd = *bndit;
    if (obnd->getBondType() != Bond::DOUBLE ||
        obnd->getStereo() <= Bond::STEREOANY ||
        obnd->getStereo() > Bond::STEREOTRANS)
      continue;

    sketcherMinimizerBondStereoInfo sinfo;
    sinfo.atom1 = ats[obnd->getStereoAtoms()[0]];
    sinfo.atom2 = ats[obnd->getStereoAtoms()[1]];
    sinfo.stereo = (obnd->getStereo() == Bond::STEREOZ ||
                    obnd->getStereo() == Bond::STEREOCIS)
                       ? sketcherMinimizerBondStereoInfo::cis
                       : sketcherMinimizerBondStereoInfo::trans;
    auto bnd = bnds[obnd->getIdx()];
    bnd->setStereoChemistry(sinfo);
    bnd->setAbsoluteStereoFromStereoInfo();
  }

  minimizer.initialize(min_mol);
  if (!params->minimizeOnly) {
    minimizer.runGenerateCoordinates();
  } else {
    CoordgenFragmenter::splitIntoFragments(min_mol);
    minimizer.m_minimizer.minimizeMolecule(min_mol);
  }
  auto conf = new Conformer(mol.getNumAtoms());
  for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    auto coords = RDGeom::Point3D(ats[i]->coordinates.x() / scaleFactor,
                                  ats[i]->coordinates.y() / scaleFactor, 0.0);
    if (params->minimizeOnly) {
      // readjust the average bond length from 1 (coordgen's target)
      // to whatever we started with
      coords *= avgBondLen;
      coords += centroid;
    }
    conf->setAtomPos(i, coords);
    // std::cerr << atom->coordinates << std::endl;
  }
  conf->set3D(false);
  mol.clearConformers();
  auto res = mol.addConformer(conf, true);
  if (!params->minimizeOnly && params->coordMap.empty() &&
      !params->templateMol) {
    // center the coordinates
    RDGeom::Transform3D tf;
    auto centroid = MolTransforms::computeCentroid(*conf);
    centroid *= -1;
    tf.SetTranslation(centroid);
    MolTransforms::transformConformer(*conf, tf);
  }
  return res;
}
}  // end of namespace CoordGen
}  // namespace RDKit
