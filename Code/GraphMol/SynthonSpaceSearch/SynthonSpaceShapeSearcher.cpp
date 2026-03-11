//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <ranges>

#include <Geometry/Transform3D.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SynthonSpaceSearch/ProgressBar.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearchHelpers.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceShapeSearcher.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

namespace RDKit::SynthonSpaceSearch {

SynthonSpaceShapeSearcher::SynthonSpaceShapeSearcher(
    const ROMol &query, const SynthonSpaceSearchParams &params,
    SynthonSpace &space)
    : SynthonSpaceSearcher(query, params, space) {
  if (space.getNumConformers() == 0) {
    throw std::runtime_error("No conformers found in SynthonSpaceSearch");
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
    const std::vector<GaussianShape::SearchShapeInput *> &fragShapes,
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
    // TODO : fix this as it's no longer true.  We only store the
    // the combination similarity.
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
          std::cout << fragShapes[fragNum]->getSmiles() << " vs "
                    << synthons[j].first << " : "
                    << synthons[j].second->getSmiles() << std::endl;
          std::cout << "sim = " << std::get<0>(sim) << " for shape "
                    << std::get<1>(sim) << std::endl;
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(j, std::get<0>(sim));
          fragMatched = true;
        }
      } else {
        unsigned int bestRefConf, bestFitConf;
        if (auto csim = fragShapes[fragNum]->bestSimilarity(
                *synthons[j].second->getShapes(), bestRefConf, bestFitConf,
                similarityCutoff);
            csim >= similarityCutoff) {
          synthonsToUse[synthonSetOrder[fragNum]][j] = true;
          fragSims[synthonSetOrder[fragNum]].emplace_back(j, csim);
          fragMatched = true;
        }
      }
    }
    if (!fragMatched) {
      std::cout << "Nothing matched " << fragShapes[fragNum]->getSmiles()
                << std::endl;
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
  std::vector<GaussianShape::SearchShapeInput *> fragShapes;
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
    for (auto &connCombPatt : connCombConnPatterns) {
      // Make sure that for this connector combination, the synthons in this
      // order have something similar.  All query fragment connectors must
      // match something in the corresponding synthon.  The synthon can
      // have unused connectors.
      bool skip = false;
      for (size_t i = 0; i < connCombPatt.size(); ++i) {
        if ((connCombPatt[i] & synthConnPatts[synthonOrder[i]]).count() <
            connCombPatt[i].count()) {
          skip = true;
          break;
        }
      }
      if (skip) {
        continue;
      }
      auto theseSynthons = getHitSynthons(
          fragShapes,
          getParams().similarityCutoff - getParams().fragSimilarityAdjuster,
          reaction, synthonOrder, *this);
      if (!theseSynthons.empty()) {
        std::unique_ptr<SynthonSpaceHitSet> hs(new SynthonSpaceShapeHitSet(
            reaction, theseSynthons, fragSet, fragShapes, synthonOrder));
        std::cout << "\nhits : " << hs->numHits << " from rxn with "
                  << reaction.getSynthons().size() << " synthons" << std::endl;
        std::cout << "input synthon order : ";
        for (auto so : synthonOrder) {
          std::cout << so << " ";
        }
        auto hst = dynamic_cast<SynthonSpaceShapeHitSet *>(hs.get());
        std::cout << "final synth order : ";
        for (size_t i = 0; i < hst->synthonSetOrder.size(); ++i) {
          std::cout << hst->synthonSetOrder[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "frags : ";
        for (const auto &f : fragSet) {
          std::cout << MolToSmiles(*f) << " : ";
        }
        std::cout << std::endl;
        if (hs->numHits) {
          results.push_back(std::move(hs));
        }
      }
    }
  }
  return results;
}

namespace {
std::unique_ptr<GaussianShape::SearchShapeInput> generateShapes(
    const ROMol &queryConfs, const ROMol &frag, double pruneThreshold) {
  // The fragSets molecules will have their atoms labelled with
  // ORIG_IDX apart from the dummy atoms, but we need coords
  // for them, too.  They are normally copied from the atom at
  // the other end of the broken bond so find that atom too.
  // Work on a copy the query.
  ROMol queryCp(queryConfs);
  std::vector<unsigned int> fragAtoms;
  fragAtoms.reserve(frag.getNumAtoms());
  boost::dynamic_bitset<> inFrag(queryCp.getNumAtoms());
  for (auto atom : frag.atoms()) {
    unsigned int origIdx;
    if (atom->getPropIfPresent<unsigned int>("ORIG_IDX", origIdx)) {
      fragAtoms.emplace_back(origIdx);
      inFrag[origIdx] = true;
    }
  }
  std::ranges::sort(fragAtoms);
  fragAtoms.erase(std::unique(fragAtoms.begin(), fragAtoms.end()),
                  fragAtoms.end());
  std::vector<std::pair<unsigned int, double>> dummyRadii;
  // std::vector<unsigned int> notColorAtoms;
  for (auto atom : frag.atoms()) {
    if (atom->getAtomicNum() == 0 && atom->getIsotope() >= 1 &&
        atom->getIsotope() <= MAX_CONNECTOR_NUM) {
      auto nbr = *frag.atomNeighbors(atom).begin();
      auto origNbr =
          queryCp.getAtomWithIdx(nbr->getProp<unsigned int>("ORIG_IDX"));
      for (auto nbrNbr : queryCp.atomNeighbors(origNbr)) {
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
  auto shapes = std::make_unique<GaussianShape::SearchShapeInput>(
      queryCp, pruneThreshold, opts);
  return shapes;
}

void generateSomeShapes(
    const std::vector<ROMol *> &fragsForShape, unsigned int beginFrag,
    unsigned int endFrag, const ROMol &queryMolHs, double pruneThreshold,
    const TimePoint *endTime,
    std::vector<std::unique_ptr<GaussianShape::SearchShapeInput>> &fragShapes) {
  if (beginFrag >= fragsForShape.size()) {
    return;
  }
  if (endFrag >= fragsForShape.size()) {
    endFrag = fragsForShape.size();
  }
  for (unsigned int fragIdx = beginFrag; fragIdx < endFrag; ++fragIdx) {
    fragShapes[fragIdx] =
        generateShapes(queryMolHs, *fragsForShape[fragIdx], pruneThreshold);
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
  std::cout << "Number of input conformations : "
            << getQuery().getNumConformers() << std::endl;
  dp_queryConfs = std::make_unique<RWMol>(getQuery(), false, -1);
  std::cout << "Number of conformers : " << dp_queryConfs->getNumConformers()
            << std::endl;
  BOOST_LOG(rdInfoLog) << "Generating query shapes for "
                       << MolToSmiles(*dp_queryConfs) << std::endl;
  std::cout << "Generating query shapes for " << MolToSmiles(*dp_queryConfs)
            << std::endl;
  std::cout << MolToCXSmiles(*dp_queryConfs) << std::endl;
  GaussianShape::ShapeInputOptions opts;
  dp_queryShapes = std::make_unique<GaussianShape::SearchShapeInput>(
      *dp_queryConfs, getParams().shapePruneThreshold, opts);
  dp_queryShapes->setActiveShape(0);

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
                                    std::ref(*dp_queryConfs),
                                    getParams().shapePruneThreshold, endTime,
                                    std::ref(d_fragShapesPool)));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    generateSomeShapes(fragsForShape, 0, fragsForShape.size(), *dp_queryConfs,
                       getParams().shapePruneThreshold, endTime,
                       d_fragShapesPool);
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
  if (!computeFragSynthonSims(endTime, minShapeSetSize)) {
    return false;
  }
#if 0
  std::cout << "Number of fragment/synthon sims : " << d_fragSynthonSims.size()
            << std::endl;
  for (const auto &fss : d_fragSynthonSims) {
    const GaussianShape::SearchShapeInput *shape =
        static_cast<const GaussianShape::SearchShapeInput *>(fss.first.first);
    const Synthon *synthon = static_cast<const Synthon *>(fss.first.second);
    std::cout << fss.first.first << " : " << shape->getSmiles()
              << " :: " << shape->getNumAtoms() << " -> " << fss.first.second
              << " : " << synthon->getSmiles()
              << " :: " << std::get<0>(fss.second) << std::endl;
  }
#endif
  return true;
}

namespace {
// Align the synth dummy and centroid with the frag equivalents and return
// the transformation that did it.
void alignDummies(const double *fragDummy, const double *fragCentroid,
                  const double *synthDummy, const double *synthCentroid,
                  RDGeom::Transform3D &xform) {
  RDGeom::Point3D fragVec =
      RDGeom::Point3D{fragDummy[0], fragDummy[1], fragDummy[2]}.directionVector(
          RDGeom::Point3D{fragCentroid[0], fragCentroid[1], fragCentroid[2]});
  RDGeom::Point3D synthVec =
      RDGeom::Point3D{synthDummy[0], synthDummy[1], synthDummy[2]}
          .directionVector(RDGeom::Point3D{synthCentroid[0], synthCentroid[1],
                                           synthCentroid[2]});
  auto axis = synthVec.crossProduct(fragVec);
  axis.normalize();
  auto cosT = fragVec.dotProduct(synthVec);
  auto sinT = sqrt(1.0 - cosT * cosT);
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

bool checkDummies(GaussianShape::ShapeInput &fragShape,
                  unsigned int fragDummyNum,
                  GaussianShape::ShapeInput &synthShape,
                  unsigned int synthDummyNum, double distThresold,
                  double angleThreshold) {
  return true;
}

// Find the best similarity of a fragment shape with a synthon shape, within
// the threshold, making sure that they start with dummy atoms aligned
// and the synthon's dummy atom doesn't move too far from the fragment's
// in the final overlay.  Assumes the fragment shape hasn't been
// normalised, so is in its original position as it was in the input
// molecule conformation.  Keeps track of the transformation of the
// best synthon shape onto the fragment for later use.
SynthonOverlay bestSimSynthonOntoFragment(
    GaussianShape::SearchShapeInput &fragShape, Synthon *synthon,
    double threshold, RDGeom::Transform3D &xform) {
  auto synthShapes = synthon->getShapes().get();
  if (fragShape.maxSimilarity(*synthShapes) < threshold) {
    return SynthonOverlay{-1.0, 0, nullptr};
  }
  GaussianShape::ShapeOverlayOptions opts;
  opts.startMode = GaussianShape::StartMode::ROTATE_0;
  opts.normalize = false;
  GaussianShape::SearchShapeInput synthShapeCp(*synthShapes);

  // Rotations of 0, 90, 180, 270 degrees.
  static const std::array<double, 4> sinT{0, sin(GaussianShape::PI / 2.0),
                                          sin(GaussianShape::PI),
                                          sin(GaussianShape::PI * 1.5)};
  static const std::array<double, 4> cosT{1.0, cos(GaussianShape::PI / 2.0),
                                          cos(GaussianShape::PI),
                                          cos(GaussianShape::PI * 1.5)};

  double bestScore = 0.0;
  unsigned int bestSynthShape = 0;
  std::shared_ptr<RDGeom::Transform3D> bestXform{new RDGeom::Transform3D()};

  std::cout << "fragShape : " << MolToCXSmiles(*fragShape.shapeToMol())
            << std::endl;
  for (unsigned int i = 0; i < fragShape.getDummyAtoms().count(); ++i) {
    unsigned int fragDummyIdx, fragDummyNbrIdx;
    fragShape.getDummyAndNbr(i, fragDummyIdx, fragDummyNbrIdx);
    auto fragDummy = fragShape.getCoords().data() + 4 * fragDummyIdx;
    auto fragDummyNbr = fragShape.getCoords().data() + 4 * fragDummyNbrIdx;
    for (unsigned int ssn = 0; ssn < synthShapes->getNumShapes(); ++ssn) {
      for (unsigned int j = 0; j < synthShapes->getDummyAtoms().count(); ++j) {
        // Align the synthon shape so that its dummy j is on top of
        // fragShape's dummy i and the vector from the dummy to its centroid
        // is aligned with the vector from the frag's dummy to centroid.
        unsigned int synthDummyIdx, synthDummyNbrIdx;
        synthShapes->getDummyAndNbr(j, synthDummyIdx, synthDummyNbrIdx);
        auto synthDummy = synthShapes->getCoords().data() + 4 * synthDummyIdx;
        auto synthDummyNbr =
            synthShapes->getCoords().data() + 4 * synthDummyNbrIdx;
        alignDummies(fragDummy, fragDummyNbr, synthDummy, synthDummyNbr, xform);
        synthShapeCp.setCoords(synthShapes->getCoords());
        synthShapeCp.transformCoords(xform);
        auto synthCpDummy = synthShapeCp.getCoords().data() + 4 * synthDummyIdx;
        auto synthCpDummyNbr =
            synthShapeCp.getCoords().data() + 4 * synthDummyNbrIdx;
#if 0
        std::cout << "synth dummy " << j << " : " << synthDummy[0] << ", "
                  << synthDummy[1] << ", " << synthDummy[2] << std::endl;
#endif
        RDGeom::Point3D rotAxis =
            RDGeom::Point3D({synthCpDummy[0], synthCpDummy[1], synthCpDummy[2]})
                .directionVector(RDGeom::Point3D(-synthCpDummyNbr[0],
                                                 -synthCpDummyNbr[1],
                                                 -synthCpDummyNbr[2]));
        RDGeom::Transform3D toOrigin;
        toOrigin.SetTranslation(RDGeom::Point3D(
            -synthCpDummy[0], -synthCpDummy[1], -synthCpDummy[2]));
        RDGeom::Transform3D fromOrigin;
        fromOrigin.SetTranslation(
            RDGeom::Point3D(synthCpDummy[0], synthCpDummy[1], synthCpDummy[2]));
        // Rotate 4 times around the axis
        for (unsigned int k = 0; k < 4; ++k) {
          // Work on a copy to leave the synthon where it started so the final
          // transform will be relevant in future.
          synthShapeCp.setCoords(synthShapes->getCoords());
          RDGeom::Transform3D rot;
          rot.SetRotation(sinT[k], cosT[k], rotAxis);
          RDGeom::Transform3D fullXform = fromOrigin * rot * toOrigin * xform;
          synthShapeCp.transformCoords(fullXform);
          RDGeom::Transform3D ovlyXform;
#if 1
          std::cout << "input singleSynthShape : "
                    << MolToCXSmiles(*synthShapeCp.shapeToMol(true))
                    << std::endl;
#endif
          auto scores = GaussianShape::AlignShape(fragShape, synthShapeCp,
                                                  &ovlyXform, opts);
#if 1
          std::cout << "Score : " << scores[0] << " " << scores[1] << " "
                    << scores[2] << std::endl;
          std::cout << "overlaid singleSynthShape : "
                    << MolToCXSmiles(*synthShapeCp.shapeToMol()) << std::endl;
          std::cout << "fragShape : " << MolToCXSmiles(*fragShape.shapeToMol())
                    << std::endl;
#endif
          if (scores[0] > bestScore) {
            bestScore = scores[0];
            // The total transform that the synthon undergoes to get this score.
            *bestXform = ovlyXform * fullXform;
            bestSynthShape = ssn;
          }
        }
      }
    }
  }
  return SynthonOverlay{bestScore, bestSynthShape, bestXform};
}

// Process a chunk of the shape->synthon similarity calculations.
void computeSomeFragSynthonSims(
    const std::vector<std::pair<GaussianShape::SearchShapeInput *, Synthon *>>
        &toDo,
    std::atomic<std::int64_t> &nextToDo, float threshold, std::mutex &mtx,
    const TimePoint *endTime, FragSynthonSims &fragSynthonSims,
    std::unique_ptr<ProgressBar> &pbar) {
  GaussianShape::SearchShapeInput *fragShape;
  Synthon *synthon;
  std::vector<float> matrix(12, 0.0);
  int numProgs = 1000;
  int numTimes = 1;
  while (true) {
    std::int64_t thisPair = ++nextToDo;
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
    fragShape = toDo[thisPair].first;
    if (!fragShape) {
      continue;
    }
    synthon = toDo[thisPair].second;
    if (synthon->getShapes()) {
      RDGeom::Transform3D xform;
      auto sim =
          bestSimSynthonOntoFragment(*fragShape, synthon, threshold, xform);
      if (std::get<0>(sim) >= threshold) {
        std::cout << fragShape->getSmiles() << " -> " << synthon->getSmiles()
                  << " : " << std::get<0>(sim)
                  << " :: " << MolToCXSmiles(*fragShape->shapeToMol()) << " :: "
                  << MolToCXSmiles(*synthon->getShapes()->shapeToMol())
                  << std::endl;
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
    const std::vector<std::pair<GaussianShape::SearchShapeInput *, Synthon *>>
        &toDo,
    float threshold, const TimePoint *endTime, FragSynthonSims &fragSynthonSims,
    std::unique_ptr<ProgressBar> &pbar, unsigned int numThreadsToUse) {
  std::atomic<std::int64_t> nextToDo = -1;
  std::vector<std::thread> threads;
  std::mutex mtx;
  if (numThreadsToUse > 1) {
    for (unsigned int i = 0u; i < numThreadsToUse; ++i) {
      threads.emplace_back(computeSomeFragSynthonSims, std::ref(toDo),
                           std::ref(nextToDo), threshold, std::ref(mtx),
                           endTime, std::ref(fragSynthonSims), std::ref(pbar));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    computeSomeFragSynthonSims(toDo, nextToDo, threshold, std::ref(mtx),
                               endTime, fragSynthonSims, pbar);
  }
}
}  // namespace

bool SynthonSpaceShapeSearcher::computeFragSynthonSims(
    const TimePoint *endTime,
    std::unordered_map<void *, unsigned int> &minShapeSetSize) {
  float threshold =
      getParams().similarityCutoff - getParams().fragSimilarityAdjuster;

  std::unique_ptr<ProgressBar> pbar;
  if (getParams().useProgressBar) {
    std::uint64_t numToDo = 0UL;
    for (size_t i = 0U; i < d_fragShapesPool.size(); ++i) {
      for (size_t j = 0U; j < getSpace().d_synthonPool.size(); ++j) {
        auto mfss = minShapeSetSize[d_fragShapesPool[i].get()];
        if (mfss <=
            getSpace().d_synthonPool[j].second->getMaxSynthonSetSize()) {
          numToDo++;
        }
      }
    }
    pbar.reset(new ProgressBar(getParams().useProgressBar, numToDo));
    std::cout << "Computing fragment/synthon shape similarities for "
              << d_fragShapesPool.size() << " fragments against "
              << getSpace().d_synthonPool.size() << " synthons with "
              << getNumThreadsToUse(getParams().numThreads) << " threads."
              << std::endl;
  }
  const auto numThreadsToUse = getNumThreadsToUse(getParams().numThreads);
  std::vector<std::pair<GaussianShape::SearchShapeInput *, Synthon *>> toDo;
  toDo.reserve(2500000);
  for (size_t i = 0U; i < d_fragShapesPool.size(); ++i) {
    for (size_t j = 0U; j < getSpace().d_synthonPool.size(); ++j) {
      auto mfss = minShapeSetSize[d_fragShapesPool[i].get()];
      // If the smallest fragment set that this fragment is in is bigger
      // than the largest SynthonSet the synthon is in, we won't ever
      // need to know the similarity between them, so skip.
      if (mfss <= getSpace().d_synthonPool[j].second->getMaxSynthonSetSize()) {
        toDo.push_back(
            std::make_pair(d_fragShapesPool[i].get(),
                           getSpace().d_synthonPool[j].second.get()));
      }
      if (toDo.size() == 2500000) {
        processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims,
                                pbar, numThreadsToUse);
        toDo.clear();
      }
    }
  }
  processShapeSynthonList(toDo, threshold, endTime, d_fragSynthonSims, pbar,
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
  double approxSim = 0.0;
  if (hasPrecomputedSims()) {
    // Calculate the overlap volumes for each fragShape->synthon using
    // the shape and colour tanimotos already calculated and use them to
    // calculate an approximation of the total similarities.
    const auto hs = dynamic_cast<const SynthonSpaceShapeHitSet *>(hitset);
    double totShapeOv = 0.0;
    double totColourOv = 0.0;
    double totSynthVol = 0.0;
    double totSynthColourVol = 0.0;
    double totFragVol = 0.0;
    double totFragColourVol = 0.0;
    // Not all fragShapes will have a matching synthon.  For example,
    // if a 2 fragment split matched 2 synthons from a 3 synthon SynthonSet.
    // All synthons of the 3rd set will have been included as potential hits.
    // Do the overlap vol calcs with the fragments to matching synthon.
    for (size_t i = 0; i < hs->fragShapes.size(); i++) {
      const auto &synth = hs->synthonsToUse[hs->synthonSetOrder[i]]
                                           [synthNums[hs->synthonSetOrder[i]]]
                                               .second;
      const auto &shapes = synth->getShapes();
      auto fragShape = hs->fragShapes[i];

      // Look up the similarity for this shape/synthon combination.  If it
      // isn't there, something has gone catastrophically wrong somewhere.
      SynthonOverlay sim;
      fragMatchedSynthon(fragShape, synth, sim);

      // Get the volume of the synthon for the shape conformer that gave
      // the similarity values.
      double synthVol = shapes->getShapeVolume(std::get<1>(sim)) -
                        shapes->getDummyVolume(std::get<1>(sim));
      double synthColourVol = shapes->getColorVolume(std::get<1>(sim));
      totSynthVol += synthVol;
      totSynthColourVol += synthColourVol;

      // There will be only shape for the fragment.
      auto fragVol = fragShape->getShapeVolume() - fragShape->getDummyVolume(0);
      auto fragColourVol = fragShape->getColorVolume();
      totFragVol += fragVol;
      totFragColourVol += fragColourVol;
      double shapeOv =
          std::get<0>(sim) * (synthVol + fragVol) / (1 + std::get<0>(sim));
      // The volumes are approximations because of the dummy volume thing, so
      // make sure the shapeOv we use isn't bigger than either of the other
      // volumes.
      shapeOv = std::min({shapeOv, synthVol, fragVol});
      totShapeOv += shapeOv;
      double colorOv = std::get<1>(sim) * (synthColourVol + fragColourVol) /
                       (1 + std::get<1>(sim));
      totColourOv += colorOv;
    }
    // Add in the volumes of un-matched Synthons
    for (size_t i = hs->fragShapes.size(); i < hs->synthonSetOrder.size();
         i++) {
      const auto &synth =
          hs->synthonsToUse[hs->synthonSetOrder[i]][synthNums[i]].second;
      const auto &shapes = synth->getShapes();
      // Use the first shape, which should have the largest volume.
      double synthVol = shapes->getShapeVolume(0) - shapes->getDummyVolume(0);
      totSynthVol += synthVol;
      double synthColourVol = shapes->getColorVolume();
      totSynthColourVol += synthColourVol;
    }
    double st = totShapeOv / (totSynthVol + totFragVol - totShapeOv);
    double ct =
        totColourOv / (totSynthColourVol + totFragColourVol - totColourOv);
    approxSim = st + ct;
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
        return false;
      }
      maxVol += shapes->getShapeVolume(0) - shapes->getDummyVolume(0);
      featureVol += shapes->getColorVolume(0);
    }
    double maxSt = std::min(maxVol, dp_queryShapes->getShapeVolume(0)) /
                   std::max(maxVol, dp_queryShapes->getShapeVolume(0));
    double maxCt = std::min(featureVol, dp_queryShapes->getColorVolume(0)) /
                   std::max(featureVol, dp_queryShapes->getColorVolume(0));
    approxSim = maxSt + maxCt;
  }
  return approxSim;
}

namespace {
// Update the hit molecule to the overlay and score represented by
// xform and sim.
void finaliseHit(
    const std::unique_ptr<RWMol> &queryConfs,
    const std::unique_ptr<RWMol> &allHitConfs,
    const std::unique_ptr<GaussianShape::SearchShapeInput> &hitShapes,
    unsigned int hitShapeNum, const RDGeom::Transform3D &xform,
    std::array<double, 3> &scores, const std::string &rxnId,
    const std::vector<const std::string *> &synthNames, ROMol &hit) {
  hit.setProp<double>("Similarity", scores[0]);
  hit.setProp<double>("ShapeTanimoto", scores[1]);
  hit.setProp<double>("ColorTanimoto", scores[2]);
  // Make a molecule of this query conformer, and add it to the hit.
  ROMol thisConf(*queryConfs, false);
  hit.setProp<std::string>("Query_CXSmiles", MolToCXSmiles(thisConf));
  const auto prodName = details::buildProductName(rxnId, synthNames);
  hit.setProp<std::string>(common_properties::_Name, prodName);

  // Copy the conformer into the hit.
  hit.clearConformers();
  auto hitConf = new Conformer(
      allHitConfs->getConformer(hitShapes->getMolConf(hitShapeNum)));
  MolTransforms::transformConformer(*hitConf, xform);
  hit.addConformer(hitConf, true);
  MolOps::assignStereochemistryFrom3D(hit);
}

// This is for when the hit was built from the hit information without
// conformation expansion, so there's less to do.  It is already overlaid.
void finaliseHit(const std::unique_ptr<RWMol> &queryConfs, ROMol &hit,
                 std::array<double, 3> &scores, const std::string &rxnId,
                 const std::vector<const std::string *> &synthNames) {
  hit.setProp<double>("Similarity", scores[0]);
  hit.setProp<double>("ShapeTanimoto", scores[1]);
  hit.setProp<double>("ColorTanimoto", scores[2]);
  // Make a molecule of this query conformer, and add it to the hit.
  ROMol thisConf(*queryConfs, false);
  hit.setProp<std::string>("Query_CXSmiles", MolToCXSmiles(thisConf));
  const auto prodName = details::buildProductName(rxnId, synthNames);
  hit.setProp<std::string>(common_properties::_Name, prodName);
}

// Transform synthon into place
std::unique_ptr<RWMol> getShapeMol(const Synthon *synthon,
                                   unsigned int shapeNum,
                                   RDGeom::Transform3D *xform) {
  // We need to build a copy of the synthon and give it the coords
  // of the associated shape that is similar to the fragment, and
  // transform it to the overlaid coords.
  const auto &synthShape = synthon->getShapes();
  GaussianShape::SearchShapeInput synthShapeCp(*synthShape);
  synthShapeCp.setActiveShape(shapeNum);
  auto shapeMol = synthShapeCp.shapeToMol(false);
  std::cout << "shape mol : " << MolToSmiles(*shapeMol, true, false, -1, false)
            << std::endl;
  std::cout << "orig mol : " << synthon->getSmiles() << std::endl;
  if (xform) {
    MolTransforms::transformConformer(shapeMol->getConformer(), *xform);
  }
  return shapeMol;
}
}  // namespace

std::unique_ptr<ROMol> SynthonSpaceShapeSearcher::buildHit(
    const SynthonSpaceHitSet *hitset, const std::vector<size_t> &synthNums,
    std::vector<const std::string *> &synthNames) const {
  std::cout << "\nBUILD HIT" << std::endl;
  // If this doesn't work, there's something so wrong an abort is necessary.
  auto hs = dynamic_cast<const SynthonSpaceShapeHitSet *>(hitset);
  PRECONDITION(hs, "Couldn't cast the hitset to shape hitset in buildHit");
  std::vector<const ROMol *> synths(synthNums.size());
  std::vector<std::shared_ptr<ROMol>> tmpSynths(synthNums.size());
  // We need to build all the synthons into a molecule, with them put into
  // the position the fragment->shape alignment gave.  Not all fragments
  // will match a synthon, such as when the reaction had 3 synthon sets
  // and there were only 2 fragments, and both those fragments had a
  // synthon that matched.  The third synthon will have been filled in from
  // all possibilities.
  boost::dynamic_bitset<> synthsMatched(synthNums.size());
  for (size_t i = 0; i < synthNums.size(); ++i) {
    const auto &synthon = hs->synthonsToUse[i][synthNums[i]].second;
    for (size_t j = 0; j < hs->fragShapes.size(); ++j) {
      const auto &shape = hs->fragShapes[j];
      // Fetch the overlay for this fragment, synthon combination
      auto it = d_fragSynthonSims.find(std::make_pair(shape, synthon));
      // This frag matched this synthon, so get the shape number and
      // transformation, and set it up.
      if (it != d_fragSynthonSims.end()) {
        unsigned int shapeNumToUse = std::get<1>(it->second);
        RDGeom::Transform3D *transToUse = std::get<2>(it->second).get();
        std::cout << "sim : " << std::get<0>(it->second)
                  << "  shapeNum : " << std::get<1>(it->second)
                  << " transform : " << std::get<2>(it->second) << std::endl;
        // We need to build a copy of the synthon and give it the coords
        // of the associated shape that is similar to the fragment, and
        // transform it to the overlaid coords.
        auto shapeMol = getShapeMol(synthon, shapeNumToUse, transToUse);
        tmpSynths[i].reset(shapeMol.release());
        std::cout << "transformed synth : " << MolToCXSmiles(*tmpSynths[i])
                  << std::endl;
        synths[i] = tmpSynths[i].get();
        synthNames[i] = &(hs->synthonsToUse[i][synthNums[i]].first);
        synthsMatched[i] = true;
      }
    }
  }
  if (synthsMatched.count() < synthNums.size()) {
    for (size_t i = 0; i < synthNums.size(); ++i) {
      if (!synthsMatched[i]) {
        std::cout << "dangling synth " << i << " : " << synthNums[i]
                  << std::endl;
        const auto &synthon = hs->synthonsToUse[i][synthNums[i]].second;
        auto shapeMol = getShapeMol(synthon, 0, nullptr);
        tmpSynths[i].reset(shapeMol.release());
        synths[i] = tmpSynths[i].get();
        synthNames[i] = &(hs->synthonsToUse[i][synthNums[i]].first);
      }
    }
  }
  return details::buildProduct(synths);
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
    auto bondlen = MolTransforms::getBondLength(conf, bond->getBeginAtomIdx(),
                                                bond->getEndAtomIdx());
    auto rad1 = PeriodicTable::getTable()->getRcovalent(
        bond->getBeginAtom()->getAtomicNum());
    auto rad2 = PeriodicTable::getTable()->getRcovalent(
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
  std::cout << "Built hit : " << MolToCXSmiles(hit) << std::endl;
  auto hitScores = GaussianShape::AlignMolecule(*dp_queryShapes, hit);
  std::cout << "hit scores : " << hitScores[0] << ", " << hitScores[1] << ", "
            << hitScores[2] << std::endl;
  std::cout << "Overlaid built hit : " << MolToCXSmiles(hit) << std::endl;
  if (checkBondLengths(hit) && !getParams().bestHit &&
      hitScores[0] >= getParams().similarityCutoff) {
    finaliseHit(dp_queryConfs, hit, hitScores, rxnId, synthNames);
    std::cout << "Initial hit is ok" << std::endl;
    return true;
  }
  std::cout << "Initial hit fails at " << hitScores[0] << std::endl;
  // If the run is multi-threaded, this will already be running
  // on the maximum number of threads, so do the embedding on
  // a single thread.
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = getParams().confRMSThreshold;
  dgParams.randomSeed = getParams().randomSeed;
  dgParams.timeout = getParams().timeOut;
  auto hitConfs =
      details::generateIsomerConformers(hit, getParams().numConformers, true,
                                        getParams().stereoEnumOpts, dgParams);
  bool foundHit = false;
  double bestSim = getParams().similarityCutoff;
  GaussianShape::ShapeInputOptions opts;
  RDGeom::Transform3D xform;
  for (auto &isomer : hitConfs) {
    // Don't prune the hit shapes, it just wastes time.
    auto hitShapes =
        std::make_unique<GaussianShape::SearchShapeInput>(*isomer, -1.0, opts);
    dp_queryShapes->setActiveShape(0);
    for (unsigned int j = 0u; j < hitShapes->getNumShapes(); ++j) {
      hitShapes->setActiveShape(j);
      auto scores =
          GaussianShape::AlignShape(*dp_queryShapes, *hitShapes, &xform);
      std::cout << "hit shape conf " << j << " of " << hitShapes->getNumShapes()
                << " : " << scores[0] << ", " << scores[1] << ", " << scores[2]
                << std::endl;
      double sim = scores[0];
      bool finalisedHit = false;
      if (sim > getBestSimilaritySoFar()) {
        finaliseHit(dp_queryConfs, isomer, hitShapes, j, xform, scores, rxnId,
                    synthNames, hit);
        finalisedHit = true;
        updateBestHitSoFar(hit, sim);
      }
      if (sim >= bestSim) {
        if (!finalisedHit) {
          finaliseHit(dp_queryConfs, isomer, hitShapes, j, xform, scores, rxnId,
                      synthNames, hit);
        }
        // If we're only interested in whether there's a shape match, and
        // not in finding the best shape, we're done.
        if (!getParams().bestHit) {
          return true;
        }
        foundHit = true;
        bestSim = scores[0];
      }
    }
  }
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
  // even though one of the other hitsets might have the smae reaction/synthon
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

  std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
      newToTry;
  newToTry.reserve(approxSims.size());
  newToTry.push_back(toTry[approxSims[0].first]);
  for (size_t i = 0; i < approxSims.size(); i++) {
    newToTry.push_back(toTry[approxSims[i].first]);
  }
  toTry = std::move(newToTry);
  details::sortAndUniquifyToTry(toTry);
  makeHitsFromToTry(toTry, endTime, results);
}

GaussianShape::SearchShapeInput *SynthonSpaceShapeSearcher::getFragShape(
    const void *frag) const {
  std::pair<const void *, ShapeSet *> tmp{frag, nullptr};
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
