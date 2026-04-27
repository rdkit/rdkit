//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <numbers>
#include <ranges>

#include <Geometry/Transform3D.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/SynthonSpaceSearch/ProgressBar.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearchHelpers.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

std::mutex myMutex;

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceShapeSearcher::SynthonSpaceShapeSearcher(
    const ROMol &query, const SynthonSpaceSearchParams &params,
    SynthonSpace *space)
    : SynthonSpaceSearcher(query, params, space) {
  if (getSpace() && getSpace()->getNumConformers() == 0) {
    throw std::runtime_error("No conformers found in SynthonSpaceSearch");
  }
  if (!getSpace()) {
    std::vector<GaussianShape::CustomFeature> allFeatures;
    buildQueryShape(query, allFeatures);
  }
}

bool SynthonSpaceShapeSearcher::fragMatchedSynthon(const void *frag,
                                                   const void *synthon,
                                                   SynthonOverlay &sim) const {
  auto it = d_fragSynthonSims.find(std::make_pair(frag, synthon));
  if (it == d_fragSynthonSims.end()) {
    sim = SynthonOverlay{-1.0, 0, nullptr};
    return false;
  }
  sim = it->second;
  return true;
}

namespace {

// Take the fragged mol ShapeSets and flag all those synthons that have a
// fragment as a similarity match.
std::vector<std::vector<size_t>> getHitSynthons(
    const std::vector<SynthonShapeInput *> &fragShapes,
    const double similarityCutoff, const SynthonSet &reaction,
    const std::vector<unsigned int> &synthonSetOrder,
    const SynthonSpaceShapeSearcher &shapeSearcher) {
  std::vector<boost::dynamic_bitset<>> synthonsToUse;
  std::vector<std::vector<size_t>> retSynthons;
  std::vector<std::vector<std::pair<size_t, double>>> fragSims(
      reaction.getSynthons().size());

  // It makes sense to match fragments against synthon sets in order of
  // smallest synthon set first because if a fragment doesn't have a match
  // in a synthon set, the whole thing's a bust.  So if fragShapes[0] is matched
  // against 1000 synthons and then fragShapes[1] is matched against 10 synthons
  // and doesn't match any of them, the first set of matches was wasted time.
  // This doesn't change the order of either, just the order of comparison.
  std::vector<std::pair<unsigned int, size_t>> fragOrders(
      synthonSetOrder.size());
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    fragOrders[i].first = i;
    fragOrders[i].second = reaction.getSynthons()[synthonSetOrder[i]].size();
  }
  std::ranges::sort(fragOrders, [](const auto &a, const auto &b) {
    return a.second < b.second;
  });
  synthonsToUse.reserve(reaction.getSynthons().size());
  for (const auto &synthonSet : reaction.getSynthons()) {
    synthonsToUse.emplace_back(synthonSet.size());
  }
  SynthonOverlay sim;
  for (size_t i = 0; i < synthonSetOrder.size(); i++) {
    const auto fragNum = fragOrders[i].first;
    const auto &synthons = reaction.getSynthons()[synthonSetOrder[fragNum]];
    // Get the smallest fragment volume.
    bool fragMatched = false;
    // Because the combination score is the sum of 2 tanimotos, it's not
    // possible to use the threshold to set upper and lower bounds on the
    // search space as is done with fingerprints and Rascal similarity.
    // We just need to plough through them in order.
    for (size_t j = 0; j < synthons.size(); j++) {
      if (!synthons[j].second->getShapes()) {
        continue;
      }
      if (shapeSearcher.hasPrecomputedSims()) {
        if (shapeSearcher.fragMatchedSynthon(fragShapes[fragNum],
                                             synthons[j].second, sim)) {
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(j, std::get<0>(sim));
          fragMatched = true;
        }
      } else {
        unsigned int bestRefConf, bestFitConf;
        RDGeom::Transform3D xform;
        if (auto csim = fragShapes[fragNum]->getShapes().bestSimilarity(
                synthons[j].second->getShapes()->getShapes(), bestRefConf,
                bestFitConf, xform, similarityCutoff);
            csim[0] >= similarityCutoff) {
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(j, csim[0]);
          fragMatched = true;
        }
      }
    }
    if (!fragMatched) {
      // No synthons matched this fragment, so the whole fragment set is a
      // bust.
      return retSynthons;
    }
  }

  // Fill in any synthons where they all didn't match because there were
  // fewer fragments than synthons.
  details::expandBitSet(synthonsToUse);
  details::bitSetsToVectors(synthonsToUse, retSynthons);

  // Now order the synthons in descending order of their similarity to
  // the corresponding fragment.
  for (size_t i = 0; i < fragShapes.size(); i++) {
    if (fragSims[i].empty()) {
      // This one will have been filled in by expandBitSet so we need to use
      // all the synthons and a dummy similarity.
      fragSims[i].resize(synthonsToUse[i].size());
      for (size_t j = 0; j < fragSims[i].size(); j++) {
        fragSims[i][j] = std::make_pair(j, 0.0);
      }
    } else {
      std::ranges::sort(
          fragSims[i].begin(), fragSims[i].end(),
          [](const auto &a, const auto &b) { return a.second > b.second; });
    }
    retSynthons[i].clear();
    retSynthons[i].reserve(fragSims[i].size());
    std::ranges::transform(fragSims[i], std::back_inserter(retSynthons[i]),
                           [](const auto &fs) { return fs.first; });
  }

  return retSynthons;
}

}  // namespace
std::vector<std::unique_ptr<SynthonSpaceHitSet>>
SynthonSpaceShapeSearcher::searchFragSet(
    const std::vector<std::shared_ptr<ROMol>> &fragSet,
    const SynthonSet &reaction) const {
  std::vector<std::unique_ptr<SynthonSpaceHitSet>> results;
  if (fragSet.size() > reaction.getSynthons().size()) {
    return results;
  }
  // Collect the ShapeSets for the fragSet
  std::vector<SynthonShapeInput *> fragShapes;
  fragShapes.reserve(fragSet.size());
  for (auto &frag : fragSet) {
    auto shape = getFragShape(frag.get());
    fragShapes.push_back(shape);
  }

  const auto connPatterns = details::getConnectorPatterns(fragSet);
  const auto synthConnPatts = reaction.getSynthonConnectorPatterns();

  // Get all the possible permutations of connector numbers compatible with
  // the number of synthon sets in this reaction.  So if the
  // fragmented molecule is C[1*].N[2*] and there are 3 synthon sets
  // we also try C[2*].N[1*], C[2*].N[3*] and C[3*].N[2*] because
  // that might be how they're labelled in the reaction database.
  const auto connCombConnPatterns =
      details::getConnectorPermutations(connPatterns, reaction.getConnectors());

  // Need to try all combinations of synthon orders.
  const auto synthonOrders =
      details::permMFromN(fragSet.size(), reaction.getSynthons().size());

  for (const auto &synthonOrder : synthonOrders) {
    auto theseSynthons = getHitSynthons(
        fragShapes,
        getParams().similarityCutoff - getParams().fragSimilarityAdjuster,
        reaction, synthonOrder, *this);
    if (!theseSynthons.empty()) {
      std::unique_ptr<SynthonSpaceHitSet> hs(new SynthonSpaceShapeHitSet(
          reaction, theseSynthons, fragSet, fragShapes, synthonOrder));
      if (hs->numHits) {
        results.push_back(std::move(hs));
      }
    }
  }
  return results;
}

namespace {
void selectFragFeatures(
    const std::vector<GaussianShape::CustomFeature> &allFeatures,
    const ROMol &queryCp, const std::vector<unsigned int> &fragAtoms,
    std::vector<GaussianShape::CustomFeature> &fragFeatures) {
  boost::dynamic_bitset<> inFrag(queryCp.getNumAtoms());
  for (const auto f : fragAtoms) {
    inFrag[f] = true;
  }

  for (const auto &feat : allFeatures) {
    bool allInFeat = true;
    // Single atom features that are dummy atoms in the query shouldn't
    // be used, but larger features, which probably just means rings,
    // should be kept.
    if (feat.atoms.size() == 1 && inFrag[feat.atoms.front()] &&
        !queryCp.getAtomWithIdx(feat.atoms.front())->getAtomicNum()) {
      allInFeat = false;
    } else {
      for (const auto a : feat.atoms) {
        if (!inFrag[a]) {
          allInFeat = false;
          break;
        }
      }
    }
    if (allInFeat) {
      fragFeatures.push_back(feat);
    }
  }
}

std::unique_ptr<SynthonShapeInput> generateShapes(
    const ROMol &queryConfs, const ROMol &frag,
    const std::vector<GaussianShape::CustomFeature> &allFeatures) {
  // The fragSets molecules will have their atoms labelled with
  // ORIG_IDX apart from the dummy atoms, but we need coords
  // for them, too.  They are normally copied from the atom at
  // the other end of the broken bond so find that atom too.
  // Work on a copy the query.
  RWMol queryCp(queryConfs);
  std::vector<unsigned int> fragAtoms;
  fragAtoms.reserve(frag.getNumAtoms());
  boost::dynamic_bitset<> inFrag(queryCp.getNumAtoms());
  for (const auto atom : frag.atoms()) {
    if (unsigned int origIdx;
        atom->getPropIfPresent<unsigned int>("ORIG_IDX", origIdx)) {
      fragAtoms.emplace_back(origIdx);
      inFrag[origIdx] = true;
    }
  }
  std::ranges::sort(fragAtoms);
  fragAtoms.erase(std::unique(fragAtoms.begin(), fragAtoms.end()),
                  fragAtoms.end());
  std::vector<std::pair<unsigned int, double>> dummyRadii;
  // std::vector<unsigned int> notColorAtoms;
  for (const auto atom : frag.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->getIsotope() >= 1 &&
        atom->getIsotope() <= MAX_CONNECTOR_NUM) {
      const auto nbr = *frag.atomNeighbors(atom).begin();
      const auto origNbr =
          queryCp.getAtomWithIdx(nbr->getProp<unsigned int>("ORIG_IDX"));
      for (const auto nbrNbr : queryCp.atomNeighbors(origNbr)) {
        if (!inFrag[nbrNbr->getIdx()]) {
          dummyRadii.emplace_back(nbrNbr->getIdx(), 2.16);
          // notColorAtoms.emplace_back(nbrNbr->getIdx());
          // Make it a dummy atom in the copy so it's correct in the Shape.
          nbrNbr->setAtomicNum(0);
        }
      }
    }
  }
  // Put the dummy atoms into the fragAtoms as we want them in the shape, too.
  // If the fragment was formed by taking a single atom out of a ring, the
  // dummy atom will appear twice in the list, so fix that.  It will mean that
  // in the shape, there will be a dummy in a ring.  A synthon would have 2, but
  // in the same place.
  std::ranges::sort(dummyRadii);
  dummyRadii.erase(std::unique(dummyRadii.begin(), dummyRadii.end()),
                   dummyRadii.end());
  std::ranges::transform(dummyRadii, std::back_inserter(fragAtoms),
                         [](const auto &p) -> unsigned int { return p.first; });

  GaussianShape::ShapeInputOptions opts;
  opts.atomSubset = fragAtoms;
  opts.atomRadii = dummyRadii;
  // Break any dummy-dummy bonds
  queryCp.beginBatchEdit();
  for (const auto bond : queryCp.bonds()) {
    if (!bond->getBeginAtom()->getAtomicNum() &&
        !bond->getEndAtom()->getAtomicNum()) {
      queryCp.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
  }
  queryCp.commitBatchEdit();
  std::vector<GaussianShape::CustomFeature> fragFeatures;
  selectFragFeatures(allFeatures, queryCp, fragAtoms, fragFeatures);
  opts.customFeatures.emplace_back(std::move(fragFeatures));
  auto shapes = std::make_unique<SynthonShapeInput>(queryCp, -1, opts);
  return shapes;
}

void generateSomeShapes(
    const std::vector<ROMol *> &fragsForShape, const unsigned int beginFrag,
    unsigned int endFrag, const ROMol &queryMol,
    const std::vector<GaussianShape::CustomFeature> &allFeatures,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<SynthonShapeInput>> &fragShapes) {
  if (beginFrag >= fragsForShape.size()) {
    return;
  }
  if (endFrag >= fragsForShape.size()) {
    endFrag = fragsForShape.size();
  }
  for (unsigned int fragIdx = beginFrag; fragIdx < endFrag; ++fragIdx) {
    fragShapes[fragIdx] =
        generateShapes(queryMol, *fragsForShape[fragIdx], allFeatures);
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    if (!(fragIdx % 100) && details::checkTimeOut(endTime)) {
      return;
    }
  }
}
}  // namespace

bool SynthonSpaceShapeSearcher::extraSearchSetup(
    std::vector<std::vector<std::shared_ptr<ROMol>>> &fragSets,
    const TimePoint *endTime) {
  // Use the given conformer unless it looks like a
  // 2D molecule.
  auto queryMol = std::unique_ptr<RWMol>(new RWMol(getQuery()));
  if (!queryMol->getNumConformers() || !queryMol->getConformer().is3D()) {
    BOOST_LOG(rdErrorLog) << "The query molecule needs a 3D conformer."
                          << std::endl;
    return false;
  }
  if (queryMol->getNumConformers() > 1) {
    BOOST_LOG(rdWarningLog)
        << "The query molecule has multiple conformers.  Just using the default."
        << std::endl;
  }
  auto queryCp = std::make_unique<RWMol>(getQuery(), false, -1);
  std::vector<GaussianShape::CustomFeature> allFeatures;
  buildQueryShape(*queryCp, allFeatures);

  // Make a map of the unique SMILES strings for the fragments, keeping
  // track of them in the vector.
  bool cancelled = false;
  auto fragSmiToFrag = details::mapFragsBySmiles(fragSets, cancelled);
  if (cancelled) {
    return false;
  }
  // Compute ShapeSets for the fragments
  d_fragShapesPool.resize(fragSmiToFrag.size());
  std::vector<ROMol *> fragsForShape;
  fragsForShape.reserve(fragSmiToFrag.size());
  std::transform(fragSmiToFrag.begin(), fragSmiToFrag.end(),
                 back_inserter(fragsForShape),
                 [](const auto &p) -> ROMol * { return p.second.front(); });

  unsigned int fragNum = 0;
  if (const auto numThreads = getNumThreadsToUse(getParams().numThreads);
      numThreads > 1) {
    const size_t eachThread = 1 + fragsForShape.size() / numThreads;
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(generateSomeShapes, std::ref(fragsForShape),
                                    start, start + eachThread,
                                    std::ref(*queryCp), allFeatures, endTime,
                                    std::ref(d_fragShapesPool)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    generateSomeShapes(fragsForShape, 0, fragsForShape.size(), *queryCp,
                       allFeatures, endTime, d_fragShapesPool);
  }
  // Keep track of the minimum fragSet size that each frag is in, which will
  // be used in computeFragSynthonSims.
  std::unordered_map<void *, unsigned int> minFragSetSize;
  for (const auto &fragSet : fragSets) {
    for (const auto &frag : fragSet) {
      if (auto it = minFragSetSize.find(frag.get());
          it != minFragSetSize.end()) {
        if (it->second < fragSet.size()) {
          it->second = fragSet.size();
        }
      } else {
        minFragSetSize[frag.get()] = fragSet.size();
      }
    }
  }
  // Use the pooled ShapeSets to populate the vectors for each fragSet.
  // Use the smallest minFragSetSize as the minShapeSetSize as we are
  // only doing 1 shape for all the identical frags.  d_fragShapes will
  // hold pointers to shapes in d_fragShapesPool, keyed on the corresponding
  // fragment.
  fragNum = 0;
  d_fragShapes.reserve(fragSmiToFrag.size());
  std::unordered_map<void *, unsigned int> minShapeSetSize;
  for (auto &[fragSmi, frags] : fragSmiToFrag) {
    minShapeSetSize[d_fragShapesPool[fragNum].get()] = 10;
    for (auto &frag : frags) {
      if (minFragSetSize[frag] <
          minShapeSetSize[d_fragShapesPool[fragNum].get()]) {
        minShapeSetSize[d_fragShapesPool[fragNum].get()] = minFragSetSize[frag];
      }
      d_fragShapes.emplace_back(frag, d_fragShapesPool[fragNum].get());
    }
    ++fragNum;
  }
  std::ranges::sort(d_fragShapes, [](const auto &p1, const auto &p2) -> bool {
    return p1.first > p2.first;
  });

  // It's likely that a fragment occurs more than once in the fragment sets,
  // so compute all the fragment->synthon similarities up front.  This makes
  // maximum use of the parallel environment, doesn't do anything that
  // wouldn't be done at some point and may prevent duplicated comparisons.
  return computeFragSynthonSims(endTime, minShapeSetSize);
}

void SynthonSpaceShapeSearcher::buildQueryShape(
    const ROMol &mol, std::vector<GaussianShape::CustomFeature> &allFeatures) {
  BOOST_LOG(rdInfoLog) << "Generating query shapes for " << MolToSmiles(mol)
                       << std::endl;
  GaussianShape::findFeatures(mol, 0, allFeatures);

  GaussianShape::ShapeInputOptions opts;
  opts.customFeatures.push_back(allFeatures);
  dp_queryShape = std::make_unique<SynthonShapeInput>(
      mol, -1, opts, getParams().shapeOverlayOptions);
}

namespace {
// Align the synth dummy and centroid with the frag equivalents and return
// the transformation that did it.
void alignDummies(const double *fragDummy, const double *fragDummyNbr,
                  const double *synthDummy, const double *synthDummyNbr,
                  RDGeom::Transform3D &xform) {
  const RDGeom::Point3D fragVec =
      RDGeom::Point3D{fragDummy[0], fragDummy[1], fragDummy[2]}.directionVector(
          RDGeom::Point3D{fragDummyNbr[0], fragDummyNbr[1], fragDummyNbr[2]});
  RDGeom::Point3D synthVec =
      RDGeom::Point3D{synthDummy[0], synthDummy[1], synthDummy[2]}
          .directionVector(RDGeom::Point3D{synthDummyNbr[0], synthDummyNbr[1],
                                           synthDummyNbr[2]});
  auto axis = synthVec.crossProduct(fragVec);
  axis.normalize();
  const auto cosT = fragVec.dotProduct(synthVec);
  const auto sinT = sqrt(1.0 - cosT * cosT);
  RDGeom::Transform3D rot;
  rot.SetRotation(cosT, sinT, axis);
  rot.TransformPoint(synthVec);

  RDGeom::Transform3D toOrigin;
  toOrigin.SetTranslation(
      RDGeom::Point3D{-synthDummy[0], -synthDummy[1], -synthDummy[2]});
  RDGeom::Transform3D fromOrigin;
  fromOrigin.SetTranslation(
      RDGeom::Point3D{fragDummy[0], fragDummy[1], fragDummy[2]});
  xform = fromOrigin * rot * toOrigin;
}

// Make sure the synthon dummy atom hasn't moved too far from the fragment
// one or the directions haven't diverged too much.
bool checkDummies(const double *fragDummyPos, const double *fragDummyNbrPos,
                  const GaussianShape::ShapeInput &synthShape,
                  const unsigned int synthDummyIdx,
                  const unsigned int synthDummyNbrIdx, double distThresholdSq,
                  const double angleThreshold) {
  const auto synthShapeDummy =
      synthShape.getCoords().data() + 3 * synthDummyIdx;
  const auto synthShapeDummyNbr =
      synthShape.getCoords().data() + 3 * synthDummyNbrIdx;

  const RDGeom::Point3D fd{fragDummyPos[0], fragDummyPos[1], fragDummyPos[2]};
  const RDGeom::Point3D sd{synthShapeDummy[0], synthShapeDummy[1],
                           synthShapeDummy[2]};
  if ((fd - sd).lengthSq() > distThresholdSq) {
    return false;
  }
  const RDGeom::Point3D fdn{fragDummyNbrPos[0], fragDummyNbrPos[1],
                            fragDummyNbrPos[2]};
  const RDGeom::Point3D sdn{synthShapeDummyNbr[0], synthShapeDummyNbr[1],
                            synthShapeDummyNbr[2]};
  const auto v1 = fd.directionVector(fdn);
  const auto v2 = sd.directionVector(sdn);
  return v1.angleTo(v2) <= angleThreshold;
}

// Find the best similarity of a fragment shape with a synthon shape, within
// the threshold, making sure that they start with dummy atoms aligned
// and the synthon's dummy atom doesn't move too far from the fragment's
// in the final overlay.  Assumes the fragment shape hasn't been
// normalised, so is in its original position as it was in the input
// molecule conformation.  Keeps track of the transformation of the
// best synthon shape onto the fragment for later use.
SynthonOverlay bestSimSynthonOntoFragment(
    const SynthonShapeInput &fragShape, const Synthon *synthon,
    double threshold,
    const GaussianShape::ShapeOverlayOptions &shapeOverlayOpts) {
  // This should only be called if the synthon has shapes.
  auto synthShapes = synthon->getShapes().get();
  if (fragShape.getShapes().maxPossibleSimilarity(synthShapes->getShapes()) <
      threshold) {
    return SynthonOverlay{-1.0, 0, nullptr};
  }
  // Work on a copy so that alpha values can be changed if need be.
  SynthonShapeInput fragShapeCp(fragShape);
  GaussianShape::ShapeOverlayOptions opts = shapeOverlayOpts;
  opts.startMode = GaussianShape::StartMode::ROTATE_0;
  opts.normalize = false;
  // Rotations of 0, 90, 180, 270 degrees.
  static const std::array<double, 4> sinT{0, sin(std::numbers::pi / 2.0),
                                          sin(std::numbers::pi),
                                          sin(std::numbers::pi * 1.5)};
  static const std::array<double, 4> cosT{1.0, cos(std::numbers::pi / 2.0),
                                          cos(std::numbers::pi),
                                          cos(std::numbers::pi * 1.5)};
  // Tolerances for how far the synthon shape moves when overlaid onto the
  // fragment.
  static constexpr double distTolSq = 4.0;  // Slightly over a bond's length
  static constexpr double angleTol =
      std::numbers::pi / 4.0;  // There's a big lever, potentially, so make
                               // the tolerance quite shallow.
  unsigned int bestSynthShape = 0;
  std::shared_ptr<RDGeom::Transform3D> bestXform{new RDGeom::Transform3D()};

  // The fragShape will only have one shape, since the query is only allowed
  // to have 1 conformation.
  // If the fragShape doesn't have any dummies, just do a normal alignment.
  if (!fragShapeCp.getNumDummyAtoms()) {
    GaussianShape::ShapeInput synthShapeCp(synthShapes->getShapes());
    unsigned int bestFragShape = 0;
    auto sim = fragShapeCp.getShapes().bestSimilarity(
        synthShapeCp, bestFragShape, bestSynthShape, *bestXform, threshold);

    return SynthonOverlay{sim[0], bestSynthShape, bestXform};
  }
  // For every dummy atom in the fragShape
  //    For every shape in the synthon:
  //       For every dummy atom in the synthon shape
  //          Align the dummy atoms and their neighbours
  //          Rotate the synthon shape by 0, 90, 190 and 270 degrees
  //          Do the alignment from that start point
  //          Keep the best score and the transformation that caused it.
  double bestScore = 0.0;
  for (const auto &[fragDummyIdx, fragDummyNbrIdx] :
       fragShapeCp.getDummyAtomsAndNbrs()) {
    auto fragDummy =
        fragShapeCp.getShapes().getCoords().data() + 3 * fragDummyIdx;
    auto fragDummyNbr =
        fragShapeCp.getShapes().getCoords().data() + 3 * fragDummyNbrIdx;
    // Set the alpha values of the other dummies -ve so they are ignored.
    // The mark missing atoms in the query so shouldn't be included in
    // any volume calculation.
    for (const auto &[otherFragDummyIdx, otherFragDummyNbr] :
         fragShapeCp.getDummyAtomsAndNbrs()) {
      if (otherFragDummyIdx != fragDummyIdx &&
          fragShapeCp.getShapes().getAlphas()[otherFragDummyIdx] > 0.0) {
        fragShapeCp.getShapes().negateAlpha(otherFragDummyIdx);
      }
    }
    for (unsigned ssn = 0; ssn < synthShapes->getShapes().getNumShapes();
         ++ssn) {
      for (const auto &[synthDummyIdx, synthDummyNbrIdx] :
           synthShapes->getDummyAtomsAndNbrs()) {
        // Fresh working copy every time.
        GaussianShape::ShapeInput singleShape(synthShapes->getShapes(), ssn);
        auto synthDummy = singleShape.getCoords().data() + 3 * synthDummyIdx;
        auto synthDummyNbr =
            singleShape.getCoords().data() + 3 * synthDummyNbrIdx;
        RDGeom::Transform3D xform;
        alignDummies(fragDummy, fragDummyNbr, synthDummy, synthDummyNbr, xform);
        singleShape.transformCoords(xform);
        auto singleShapeDummy =
            singleShape.getCoords().data() + 3 * synthDummyIdx;
        auto singleShapeDummyNbr =
            singleShape.getCoords().data() + 3 * synthDummyNbrIdx;
        RDGeom::Point3D rotAxis =
            RDGeom::Point3D(
                {singleShapeDummy[0], singleShapeDummy[1], singleShapeDummy[2]})
                .directionVector(RDGeom::Point3D(singleShapeDummyNbr[0],
                                                 singleShapeDummyNbr[1],
                                                 singleShapeDummyNbr[2]));
        // Get the new dummy->neighbour axis to rotate about.
        RDGeom::Transform3D toOrigin;
        toOrigin.SetTranslation(RDGeom::Point3D(
            -singleShapeDummy[0], -singleShapeDummy[1], -singleShapeDummy[2]));
        RDGeom::Transform3D fromOrigin;
        fromOrigin.SetTranslation(RDGeom::Point3D(
            singleShapeDummy[0], singleShapeDummy[1], singleShapeDummy[2]));
        // Rotate 4 times around the axis
        for (unsigned int k = 0; k < 4; ++k) {
          // Work on a copy so the transformation will be relevant at the end.
          GaussianShape::ShapeInput singleShapeCp = singleShape;
          RDGeom::Transform3D rot;
          rot.SetRotation(sinT[k], cosT[k], rotAxis);
          RDGeom::Transform3D fullXform = fromOrigin * rot * toOrigin;
          singleShapeCp.transformCoords(fullXform);
          RDGeom::Transform3D ovlyXform;
          auto scores = GaussianShape::AlignShape(
              fragShapeCp.getShapes(), singleShapeCp, &ovlyXform, opts);
          if (scores[0] > bestScore) {
            // Reject anything where the synthon dummies have gone too far.
            if (checkDummies(fragDummy, fragDummyNbr, singleShapeCp,
                             synthDummyIdx, synthDummyNbrIdx, distTolSq,
                             angleTol)) {
              bestScore = scores[0];
              // The total transform that the synthon undergoes to get this
              // score.
              *bestXform = ovlyXform * fullXform * xform;
              bestSynthShape = ssn;
            }
          }
        }
      }
    }
    // Set the alpha values of the other dummies back again.
    for (const auto &[otherFragDummyIdx, otherFragDummyNbr] :
         fragShapeCp.getDummyAtomsAndNbrs()) {
      if (otherFragDummyIdx != fragDummyIdx &&
          fragShapeCp.getShapes().getAlphas()[otherFragDummyIdx] < 0.0) {
        fragShapeCp.getShapes().negateAlpha(otherFragDummyIdx);
      }
    }
  }
  return SynthonOverlay{bestScore, bestSynthShape, bestXform};
}

// Process a chunk of the shape->synthon similarity calculations.
void computeSomeFragSynthonSims(
    const std::vector<std::pair<SynthonShapeInput *, Synthon *>> &toDo,
    std::atomic<std::int64_t> &nextToDo, float threshold, std::mutex &mtx,
    const TimePoint *endTime, FragSynthonSims &fragSynthonSims,
    const GaussianShape::ShapeOverlayOptions &shapeOverlayOpts,
    std::unique_ptr<ProgressBar> &pbar) {
  std::vector<float> matrix(12, 0.0);
  int numProgs = 1000;
  int numTimes = 1;
  while (true) {
    const std::int64_t thisPair = ++nextToDo;
    if (thisPair >= static_cast<std::int64_t>(toDo.size())) {
      return;
    }
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    --numTimes;
    if (!numTimes) {
      numTimes = 100;
      if (details::checkTimeOut(endTime)) {
        return;
      }
    }
    SynthonShapeInput *fragShape = toDo[thisPair].first;
    if (!fragShape) {
      continue;
    }
    Synthon *synthon = toDo[thisPair].second;
    if (synthon->getShapes()) {
      auto sim = bestSimSynthonOntoFragment(*fragShape, synthon, threshold,
                                            shapeOverlayOpts);
      if (std::get<0>(sim) >= threshold) {
        std::pair<const void *, const void *> p{fragShape, synthon};
        std::unique_lock lock1{mtx};
        fragSynthonSims.insert(std::make_pair(p, sim));
      }
      if (pbar) {
        --numProgs;
        if (!numProgs) {
          numProgs = 1000;
          pbar->increment(1000);
        }
      }
    }
  }
}

// Calculate shape->synthon similarities for this block of the full search
// space, in parallel if required.
void processShapeSynthonList(
    const std::vector<std::pair<SynthonShapeInput *, Synthon *>> &toDo,
    float threshold, const TimePoint *endTime, FragSynthonSims &fragSynthonSims,
    const GaussianShape::ShapeOverlayOptions &shapeOverlayOpts,
    std::unique_ptr<ProgressBar> &pbar, const unsigned int numThreadsToUse) {
  std::atomic<std::int64_t> nextToDo = -1;
  std::mutex mtx;
  if (numThreadsToUse > 1) {
    std::vector<std::thread> threads;
    threads.reserve(numThreadsToUse);
    for (unsigned int i = 0u; i < numThreadsToUse; ++i) {
      threads.emplace_back(computeSomeFragSynthonSims, std::ref(toDo),
                           std::ref(nextToDo), threshold, std::ref(mtx),
                           endTime, std::ref(fragSynthonSims),
                           std::ref(shapeOverlayOpts), std::ref(pbar));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    computeSomeFragSynthonSims(toDo, nextToDo, threshold, std::ref(mtx),
                               endTime, fragSynthonSims, shapeOverlayOpts,
                               pbar);
  }
}
}  // namespace

bool SynthonSpaceShapeSearcher::computeFragSynthonSims(
    const TimePoint *endTime,
    std::unordered_map<void *, unsigned int> &minFragSetSize) {
  const float threshold =
      getParams().similarityCutoff - getParams().fragSimilarityAdjuster;

  std::unique_ptr<ProgressBar> pbar;
  if (getParams().useProgressBar) {
    std::uint64_t numToDo = 0UL;
    for (const auto &fragShape : d_fragShapesPool) {
      for (const auto &synthon : getSpace()->d_synthonPool) {
        const auto mfss = minFragSetSize[fragShape.get()];
        if (mfss <= synthon.second->getMaxSynthonSetSize()) {
          numToDo++;
        }
      }
    }
    pbar.reset(new ProgressBar(getParams().useProgressBar, numToDo));
    std::cout << "Computing fragment/synthon shape similarities for "
              << d_fragShapesPool.size() << " fragments against "
              << getSpace()->d_synthonPool.size() << " synthons with "
              << getNumThreadsToUse(getParams().numThreads) << " threads."
              << std::endl;
  }
  const auto numThreadsToUse = getNumThreadsToUse(getParams().numThreads);
  std::vector<std::pair<SynthonShapeInput *, Synthon *>> toDo;
  toDo.reserve(2500000);
  for (const auto &fragShape : d_fragShapesPool) {
    for (const auto &synthon : getSpace()->d_synthonPool) {
      const auto mfss = minFragSetSize[fragShape.get()];
      // If the smallest fragment set that this fragment is in is bigger
      // than the largest SynthonSet the synthon is in, we won't ever
      // need to know the similarity between them, so skip.
      if (mfss <= synthon.second->getMaxSynthonSetSize()) {
        toDo.push_back(std::make_pair(fragShape.get(), synthon.second.get()));
      }
      if (toDo.size() == 2500000) {
        processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims,
                                getParams().shapeOverlayOptions, pbar,
                                numThreadsToUse);
        toDo.clear();
      }
    }
  }
  processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims,
                          getParams().shapeOverlayOptions, pbar,
                          numThreadsToUse);
  return !ControlCHandler::getGotSignal();
}  // namespace

bool SynthonSpaceShapeSearcher::quickVerify(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  // The approximate similarity should already have been checked in
  // processToTrySet so no need to do it again.
  if (!SynthonSpaceSearcher::quickVerify(hitset, synthNums)) {
    return false;
  }
  return true;
}

double SynthonSpaceShapeSearcher::approxSimilarity(
    const SynthonSpaceHitSet *hitset,
    const std::vector<size_t> &synthNums) const {
  // If this doesn't work, there's something so wrong an abort is necessary.
  const auto hs = dynamic_cast<const SynthonSpaceShapeHitSet *>(hitset);
  PRECONDITION(hs, "Couldn't cast the hitset to shape hitset in buildHit");

  const auto &sso = hs->synthonSetOrder;
  double approxSim = 0.0;
  if (hasPrecomputedSims()) {
    // Take the approximate similarity as the mean of the individual
    // frag->synthon similarities weighted by the synthon volumes.
    double totVol = 0.0;
    std::vector<double> synthVols(synthNums.size(), 0.0);
    std::vector<double> scores(synthNums.size(), 0.0);
    auto &m = getParams().shapeOverlayOptions.optParam;

    // Not all fragShapes will have a matching synthon.  For example,
    // if a 2 fragment split matched 2 synthons from a 3 synthon SynthonSet.
    // All synthons of the 3rd set will have been included as potential hits.
    // Do the overlap vol calcs with the fragments to matching synthon.
    boost::dynamic_bitset<> synthsMatched(synthNums.size());
    for (size_t i = 0; i < hs->fragShapes.size(); ++i) {
      const auto &synthon = hs->synthonsToUse[sso[i]][synthNums[sso[i]]].second;
      const auto &shapes = synthon->getShapes();
      if (!shapes) {
        // If one of the synthons doesn't have a shape, the hit is clearly
        // a bust.
        return 0.0;
      }
      const auto &fragShape = hs->fragShapes[i];
      // Look up the similarity for this shape/synthon combination.  If it
      // isn't there, something has gone catastrophically wrong somewhere.
      SynthonOverlay sim;
      if (fragMatchedSynthon(fragShape, synthon, sim)) {
        // Get the volume of the synthon for the shape conformer that gave
        // the similarity values.
        const double synthVol =
            shapes->getShapes().getShapeVolume(std::get<1>(sim)) -
            shapes->getDummyVolume(std::get<1>(sim));
        const double synthColourVol =
            shapes->getShapes().getColorVolume(std::get<1>(sim));

        synthVols[sso[i]] = (1.0 - m) * synthVol + m * synthColourVol;
        scores[sso[i]] = std::get<0>(sim);
        totVol += synthVols[sso[i]];
        synthsMatched[sso[i]] = true;
      } else {
        return 0.0;
      }
    }
    if (synthsMatched.count() < synthNums.size()) {
      for (size_t i = 0; i < synthNums.size(); ++i) {
        if (!synthsMatched[i]) {
          const auto &synthon = hs->synthonsToUse[i][synthNums[i]].second;
          // Use the first shape, which should have the largest volume.
          const auto &shapes = synthon->getShapes();
          if (!shapes) {
            return 0.0;
          }
          const double synthVol =
              shapes->getShapes().getShapeVolume(0) - shapes->getDummyVolume(0);
          const double synthColourVol = shapes->getShapes().getColorVolume(0);
          synthVols[i] = (1.0 - m) * synthVol + m * synthColourVol;
          totVol += synthVols[i];
        }
      }
    }
    // scores[i] will be zero for un-matched synthons.
    for (unsigned int i = 0; i < synthNums.size(); i++) {
      approxSim += scores[i] * synthVols[i] / totVol;
    }
  } else {
    // This is the best we can do without pre-computed similarities
    double maxVol = 0.0;
    double featureVol = 0.0;
    // The synthon shapes are sorted in descending order of sov + sof.
    // Assume therefore that the maximum volume of the hit is the sum
    // of the sov's of the first shape in each synthon, minus the volume
    // of their dummy atoms.
    for (size_t i = 0; i < synthNums.size(); i++) {
      const auto &synth = hitset->synthonsToUse[i][synthNums[i]].second;
      const auto &shapes = synth->getShapes();
      if (!shapes) {
        return 0.0;
      }
      maxVol +=
          shapes->getShapes().getShapeVolume(0) - shapes->getDummyVolume(0);
      featureVol += shapes->getShapes().getColorVolume(0);
    }
    approxSim = GaussianShape::maxScore(
        maxVol, dp_queryShape->getShapes().getShapeVolume(0), featureVol,
        dp_queryShape->getShapes().getColorVolume(0),
        getParams().shapeOverlayOptions);
  }
  return approxSim;
}

namespace {
double calcExcludedVolume(const ROMol &mol,
                          const GaussianShape::ShapeInput &excVol,
                          const GaussianShape::ShapeOverlayOptions &ovlyOpts) {
  GaussianShape::ShapeInput hitShape(
      mol, -1, GaussianShape::ShapeInputOptions(), ovlyOpts);
  std::array<double, 2> excVols;
  GaussianShape::ScoreShape(hitShape, excVol, ovlyOpts, &excVols);
  return excVols[0] + excVols[1];
}

unsigned int calcNumClashes(const ROMol &mol,
                            const GaussianShape::ShapeInput &excVol) {
  static constexpr double doubleCsq =
      4 * GaussianShape::CARBON_RAD * GaussianShape::CARBON_RAD;
  const auto shpCoords = excVol.getCoords();
  boost::dynamic_bitset<> clashAtoms(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    auto aPos = mol.getConformer().getAtomPos(atom->getIdx());
    for (unsigned int i = 0; i < shpCoords.size(); i += 3) {
      const RDGeom::Point3D sPos{shpCoords[i], shpCoords[i + 1],
                                 shpCoords[i + 2]};
      if ((aPos - sPos).lengthSq() < doubleCsq) {
        clashAtoms[atom->getIdx()] = true;
        break;
      }
    }
  }
  return clashAtoms.count();
}

bool checkExcludedVols(const ROMol &mol, const SynthonSpaceSearchParams &params,
                       double &excludedVol, double &meanExcludedVol) {
  if (params.excludedVolume) {
    excludedVol = calcExcludedVolume(mol, *params.excludedVolume,
                                     params.shapeOverlayOptions);
    auto numClashes = calcNumClashes(mol, *params.excludedVolume);
    meanExcludedVol = excludedVol / static_cast<double>(numClashes);
    if (params.maxExcludedVolume > -0.5 &&
        excludedVol > params.maxExcludedVolume) {
      return false;
    }
    if (params.maxMeanExcludedVolume > -0.5) {
      if (meanExcludedVol > params.maxMeanExcludedVolume) {
        return false;
      }
    }
  }
  return true;
}

// Assuming the hitShape is aligned onto the query, extract that conformation
// into hit and add the scores to it as properties..  Return false if it
// fails one of the excluded volume checks, true otherwise.
bool finaliseHit(GaussianShape::ShapeInput &hitShapes,
                 const unsigned int hitShapeNum,
                 const std::array<double, 3> &scores, const std::string &rxnId,
                 const std::vector<const std::string *> &synthNames, ROMol &hit,
                 const SynthonSpaceSearchParams &params, double &excludedVol,
                 double &meanExcludedVol) {
  // Copy the conformer into the hit.
  hitShapes.setActiveShape(hitShapeNum);
  hit = static_cast<ROMol>(*hitShapes.shapeToMol(false, true));
  if (!checkExcludedVols(hit, params, excludedVol, meanExcludedVol)) {
    return false;
  }
  hit.setProp<double>("Similarity", scores[0]);
  hit.setProp<double>("ShapeScore", scores[1]);
  hit.setProp<double>("ColorScore", scores[2]);
  const auto prodName = details::buildProductName(rxnId, synthNames);
  hit.setProp<std::string>(common_properties::_Name, prodName);
  MolOps::assignStereochemistryFrom3D(hit);
  if (excludedVol > -0.5) {
    hit.setProp<double>("ExcludedVolume", excludedVol);
  }
  if (meanExcludedVol > -0.5) {
    hit.setProp<double>("MeanExcludedVolume", meanExcludedVol);
  }
  return true;
}

// This is for when the hit was built from the hit information without
// conformation expansion, so there's less to do.  It is already overlaid,
// so just add the scores as properties.
void finaliseHit(const ROMol &hit, const std::array<double, 3> &scores,
                 const std::string &rxnId,
                 const std::vector<const std::string *> &synthNames,
                 const double excludedVol, const double meanExcludedVol) {
  hit.setProp<double>("Similarity", scores[0]);
  hit.setProp<double>("ShapeScore", scores[1]);
  hit.setProp<double>("ColorScore", scores[2]);
  const auto prodName = details::buildProductName(rxnId, synthNames);
  hit.setProp<std::string>(common_properties::_Name, prodName);
  if (excludedVol > -0.5) {
    hit.setProp<double>("ExcludedVolume", excludedVol);
  }
  if (meanExcludedVol > -0.5) {
    hit.setProp<double>("MeanExcludedVolume", meanExcludedVol);
  }
}

// Transform synthon into place
std::unique_ptr<RWMol> getShapeMol(const Synthon *synthon,
                                   const unsigned int shapeNum,
                                   const RDGeom::Transform3D *xform) {
  // Make a molecule from the synthon's shape, and
  // transform it to the overlaid coords.
  const auto &synthShape = synthon->getShapes();
  synthShape->getShapes().setActiveShape(shapeNum);
  auto shapeMol = synthon->getShapes()->getShapes().shapeToMol(false, true);
  if (xform) {
    MolTransforms::transformConformer(shapeMol->getConformer(), *xform);
  }
  return shapeMol;
}
}  // namespace

std::unique_ptr<ROMol> SynthonSpaceShapeSearcher::buildHit(
    const SynthonSpaceHitSet *hitset, const std::vector<size_t> &synthNums,
    std::vector<const std::string *> &synthNames) const {
  // If this doesn't work, there's something so wrong an abort is necessary.
  const auto hs = dynamic_cast<const SynthonSpaceShapeHitSet *>(hitset);
  PRECONDITION(hs, "Couldn't cast the hitset to shape hitset in buildHit");
  // shapeSynths is a molecule produced by shapeFromMol, with coords.
  std::vector<const ROMol *> shapeSynths(synthNums.size());
  // plainSynths is the original synthon input mol, to be used if the
  // shapeSynths molecules can't be zipped.  The commonest reason for that
  // is the shape was built from a sampleMol that didn't have exactly the
  // same chemistry.
  std::vector<const ROMol *> plainSynths(synthNums.size());
  std::vector<std::shared_ptr<ROMol>> tmpSynths(synthNums.size());
  // We need to build all the synthons into a molecule, with them put into
  // the position the fragment->shape alignment gave.  Not all fragments
  // will match a synthon, such as when the reaction had 3 synthon sets
  // and there were only 2 fragments, and both those fragments had a
  // synthon that matched.  The third synthon will have been filled in from
  // all possibilities.
  boost::dynamic_bitset<> synthsMatched(synthNums.size());
  const auto &sso = hs->synthonSetOrder;
  for (size_t i = 0; i < hs->fragShapes.size(); ++i) {
    const auto &synthon = hs->synthonsToUse[sso[i]][synthNums[sso[i]]].second;
    const auto &shapes = synthon->getShapes();
    if (!shapes) {
      // If one of the synthons doesn't have a shape, the hit is clearly
      // a bust.
      return std::unique_ptr<ROMol>();
    }
    const auto &fragShape = hs->fragShapes[i];
    // Look up the similarity for this shape/synthon combination.  If it
    // isn't there, something has gone catastrophically wrong somewhere.
    SynthonOverlay sim;
    if (fragMatchedSynthon(fragShape, synthon, sim)) {
      const unsigned int shapeNumToUse = std::get<1>(sim);
      const RDGeom::Transform3D *transToUse = std::get<2>(sim).get();
      // We need to build a copy of the synthon and give it the coords
      // of the associated shape that is similar to the fragment, and
      // transform it to the overlaid coords.
      auto shapeMol = getShapeMol(synthon, shapeNumToUse, transToUse);
      tmpSynths[sso[i]].reset(shapeMol.release());
      shapeSynths[sso[i]] = tmpSynths[sso[i]].get();
      plainSynths[sso[i]] = synthon->getOrigMol().get();
      synthNames[sso[i]] = &hs->synthonsToUse[sso[i]][synthNums[sso[i]]].first;
      synthsMatched[sso[i]] = true;
    }
  }
  if (synthsMatched.count() < synthNums.size()) {
    for (size_t i = 0; i < synthNums.size(); ++i) {
      if (!synthsMatched[i]) {
        const auto &synthon = hs->synthonsToUse[i][synthNums[i]].second;
        auto shapeMol = getShapeMol(synthon, 0, nullptr);
        tmpSynths[i].reset(shapeMol.release());
        shapeSynths[i] = tmpSynths[i].get();
        synthNames[i] = &hs->synthonsToUse[i][synthNums[i]].first;
        plainSynths[i] = synthon->getOrigMol().get();
      }
    }
  }
  auto prod = details::buildProduct(shapeSynths);
  if (!prod) {
    // Try it with the plain synthons, which will obviously need coordinates
    // creating later on.
    prod = details::buildProduct(plainSynths);
    if (!prod) {
      std::cout << "Failed building " << std::endl;
      for (const auto &s : synthNames) {
        std::cout << *s << " ";
      }
      std::cout << " for reaction " << hs->d_reaction->getId() << std::endl;
    }
  }
  return prod;
}

namespace {
bool checkBondLengths(const ROMol &mol) {
  // DetermineBonds::connectivityVdw uses a covalent factor of 1.3.
  static constexpr double radFactor = 1.3;
  const auto conf = mol.getConformer();
  for (const auto bond : mol.bonds()) {
    if (!bond->getBeginAtom()->getAtomicNum() ||
        !bond->getEndAtom()->getAtomicNum()) {
      continue;
    }
    const auto bondlen = MolTransforms::getBondLength(
        conf, bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    const auto rad1 = PeriodicTable::getTable()->getRcovalent(
        bond->getBeginAtom()->getAtomicNum());
    const auto rad2 = PeriodicTable::getTable()->getRcovalent(
        bond->getEndAtom()->getAtomicNum());
    if (bondlen > radFactor * (rad1 + rad2)) {
      return false;
    }
  }
  return true;
}
}  // namespace

bool SynthonSpaceShapeSearcher::verifyHit(
    ROMol &hit, const std::string &rxnId,
    const std::vector<const std::string *> &synthNames) {
  // std::cout << "Verifying hit " << MolToCXSmiles(hit) << std::endl;
  std::array<double, 3> initScores{-1.0, -1.0, -1.0};
  if (hit.getNumConformers()) {
    initScores = GaussianShape::AlignMolecule(
        dp_queryShape->getShapes(), hit, GaussianShape::ShapeInputOptions(),
        nullptr, getParams().shapeOverlayOptions);
    double excludedVol{-1.0}, meanExcludedVol{-1.0};
    if (checkExcludedVols(hit, getParams(), excludedVol, meanExcludedVol)) {
      // If the excluded vol is also ok, take this hit.  Otherwise, pass it
      // through to conformation expansion to see if it wins there.
      finaliseHit(hit, initScores, rxnId, synthNames, excludedVol,
                  meanExcludedVol);
      if (!getParams().bestHit && checkBondLengths(hit) &&
          initScores[0] >= getParams().similarityCutoff) {
        return true;
      }
    }
  }
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.pruneRmsThresh = getParams().confRMSThreshold;
  dgParams.randomSeed = getParams().randomSeed;
  dgParams.timeout = getParams().timeOut;
  UserConfGenerator userConf = getParams().userConformerGenerator;
  auto hitConfs = details::generateIsomerConformers(
      hit, getParams().numConformers, true, getParams().stereoEnumOpts,
      dgParams, userConf);
  bool foundHit = false;
  // Prefer a hit from a generated conf rather than the initial hit
  // which might have iffy bond lengths because it's the output from molzip.
  double bestSim = getParams().similarityCutoff;
  GaussianShape::ShapeInputOptions opts;
  GaussianShape::ShapeOverlayOptions overlayOptions =
      getParams().shapeOverlayOptions;
  // Don't prune the hit shapes, it just wastes time.
  opts.shapePruneThreshold = -1.0;
  RDGeom::Transform3D xform;
  for (auto &isomer : hitConfs) {
    GaussianShape::ShapeInput isomerShapes(*isomer, -1, opts, overlayOptions);
    RDGeom::Transform3D bestXform;
    // Do all the shapes one by one rather than via bestSimilarity because
    // excluded volumes might be too high for the best scoring overlay, but
    // there might be an acceptable one with a lower score.
    for (unsigned int shapeNum = 0; shapeNum < isomerShapes.getNumShapes();
         shapeNum++) {
      double excludedVol{-1.0}, meanExcludedVol{-1.0};
      bool finalisedHit = false;
      isomerShapes.setActiveShape(shapeNum);
      auto thisScore = GaussianShape::AlignShape(
          dp_queryShape->getShapes(), isomerShapes, &bestXform, overlayOptions);
      if (thisScore[0] > getBestSimilaritySoFar()) {
        if (!finaliseHit(isomerShapes, shapeNum, thisScore, rxnId, synthNames,
                         hit, getParams(), excludedVol, meanExcludedVol)) {
          // It failed excluded vols so don't use it.
          continue;
        }
        finalisedHit = true;
        std::unique_lock lock1{myMutex};
        updateBestHitSoFar(hit, thisScore[0]);
      }
      if (thisScore[0] >= bestSim &&
          thisScore[0] > getParams().similarityCutoff) {
        if (!finalisedHit) {
          if (!finaliseHit(isomerShapes, shapeNum, thisScore, rxnId, synthNames,
                           hit, getParams(), excludedVol, meanExcludedVol)) {
            continue;
          }
        }
        // If we're only interested in whether there's a shape match, and
        // not in finding the best shape, we're done.
        if (!getParams().bestHit) {
          return true;
        }
        foundHit = true;
        bestSim = thisScore[0];
      }
      if (!foundHit && initScores[0] > getParams().similarityCutoff) {
        // Stick with what we found at the start, which should still be in hit.
        foundHit = true;
      }
    }
  }
  // std::cout << "Hit : " << MolToCXSmiles(hit) << std::endl;
  // std::cout << "Hit similarity : " << hit.getProp<double>("Similarity") <<
  // " : "
  //           << hit.getProp<std::string>(common_properties::_Name) <<
  //           std::endl;
  return foundHit;
}

void SynthonSpaceShapeSearcher::processToTrySet(
    std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
        &toTry,
    const TimePoint *endTime, std::vector<std::unique_ptr<ROMol>> &results) {
  // This needs to be done a bit differently because the approximate
  // similarities depend on the fragments as well as the synthons selected.
  // For the other searches, such as fingerprints, the synthons are compared
  // with the full query so it doesn't matter which entry in the hit set
  // is used.  If we uniquify on reaction/synthons straight away, we
  // may end up with a reaction/synthon pair that doesn't hit the threshold
  // even though one of the other hitsets might have the same reaction/synthon
  // pair coupled with a different fragSet that does go over the threshold.
  // First calculate the approximate similarities and discard any that
  // don't meet the threshold.
  std::vector<std::pair<size_t, double>> approxSims;
  approxSims.reserve(toTry.size());
  for (size_t i = 0U; i < toTry.size(); ++i) {
    const auto &tt = toTry[i];
    auto sim = approxSimilarity(tt.first, tt.second);
    if (sim >=
        getParams().similarityCutoff - getParams().approxSimilarityAdjuster) {
      approxSims.push_back(std::make_pair(i, sim));
    }
  }
  if (approxSims.empty()) {
    return;
  }
  std::sort(approxSims.begin(), approxSims.end(),
            [](const auto &p1, const auto &p2) -> bool {
              return p1.second > p2.second;
            });
  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
      newToTry;
  newToTry.reserve(approxSims.size());

  for (const auto &as : approxSims) {
    newToTry.push_back(toTry[as.first]);
  }
  toTry = std::move(newToTry);
  details::sortAndUniquifyToTry(toTry);
  makeHitsFromToTry(toTry, endTime, results);
}

SynthonShapeInput *SynthonSpaceShapeSearcher::getFragShape(
    const void *frag) const {
  const std::pair<const void *, ShapeSet *> tmp{frag, nullptr};
  const auto it = std::ranges::lower_bound(
      d_fragShapes, tmp, [](const auto &p1, const auto &p2) -> bool {
        return p1.first > p2.first;
      });
  if (it != d_fragShapes.end()) {
    return it->second;
  }
  return nullptr;
}

}  // namespace RDKit::SynthonSpaceSearch
