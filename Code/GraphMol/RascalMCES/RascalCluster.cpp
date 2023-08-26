//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains an implementation of the clustering algorithm
// described in
// 'A Line Graph Algorithm for Clustering Chemical Structures Based
// on Common Substructural Cores', JW Raymond, PW Willett.
// https://match.pmf.kg.ac.rs/electronic_versions/Match48/match48_197-207.pdf
// https://eprints.whiterose.ac.uk/77598/
// It uses the RASCAL MCES algorithm to perform a fuzzy clustering
// of a set of molecules.

#include <algorithm>
#include <iterator>
#include <list>
#include <thread>
#include <vector>

#include <RDGeneral/RDThreads.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalDetails.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalResult.h>

namespace RDKit {
namespace RascalMCES {
namespace details {
ClusNode calcMolMolSimilarity(
    const std::tuple<
        size_t, size_t, const std::vector<std::shared_ptr<ROMol>> *,
        const RascalOptions *, const RascalClusterOptions *> &toDo) {
  auto i = std::get<0>(toDo);
  auto j = std::get<1>(toDo);
  auto mols = std::get<2>(toDo);
  auto opts = std::get<3>(toDo);
  auto clusOpts = std::get<4>(toDo);
  auto res = rascalMCES(*(*mols)[i], *(*mols)[j], *opts);
  ClusNode cn;
  cn.d_mol1Num = i;
  cn.d_mol2Num = j;
  if (res.empty()) {
    // tier1Sim and tier2Sim were above the threshold, but no MCES
    // was found.
    cn.d_sim = 0.0;
  } else {
    if (res.front().getBondMatches().empty()) {
      cn.d_sim = 0.0;
    } else {
      res.front().trimSmallFrags();
      res.front().largestFragsOnly(clusOpts->maxNumFrags);
      cn.d_sim = res.front().getSimilarity();
      if (cn.d_sim >= opts->similarityThreshold) {
        cn.d_res = std::shared_ptr<RascalResult>(new RascalResult(res.front()));
      }
    }
  }
  return cn;
}

std::vector<std::vector<ClusNode>> buildProximityGraph(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts) {
  if (mols.size() < 2) {
    return std::vector<std::vector<ClusNode>>();
  }
  std::vector<std::vector<ClusNode>> proxGraph =
      std::vector<std::vector<ClusNode>>(
          mols.size(), std::vector<ClusNode>(mols.size(), ClusNode()));
  std::vector<
      std::tuple<size_t, size_t, const std::vector<std::shared_ptr<ROMol>> *,
                 const RascalOptions *, const RascalClusterOptions *>>
      toDo;

  RascalOptions opts;
  opts.similarityThreshold = clusOpts.similarityCutoff;
  for (size_t i = 0; i < mols.size() - 1; ++i) {
    for (size_t j = i + 1; j < mols.size(); ++j) {
      toDo.push_back({i, j, &mols, &opts, &clusOpts});
    }
  }

  auto buildProxGraphPart =
      [](const std::vector<std::tuple<
             size_t, size_t, const std::vector<std::shared_ptr<ROMol>> *,
             const RascalOptions *, const RascalClusterOptions *>> &toDo,
         std::vector<ClusNode> &molSims, size_t start, size_t finish) -> void {
    if (start > toDo.size()) {
      return;
    }
    if (finish > toDo.size()) {
      finish = toDo.size();
    }
    std::transform(toDo.begin() + start, toDo.begin() + finish,
                   molSims.begin() + start, calcMolMolSimilarity);
  };

  std::vector<ClusNode> molSims(toDo.size());
#if RDK_BUILD_THREADSAFE_SSS
  auto numThreads = getNumThreadsToUse(clusOpts.numThreads);
  if (numThreads > 1) {
    size_t eachThread = 1 + (toDo.size() / numThreads);
    size_t start = 0;
    std::vector<std::thread> threads;
    for (unsigned int i = 0U; i < numThreads; ++i, start += eachThread) {
      threads.push_back(std::thread(buildProxGraphPart, std::ref(toDo),
                                    std::ref(molSims), start,
                                    start + eachThread));
    }
    for (auto &t : threads) {
      t.join();
    }
  } else {
    std::transform(toDo.begin(), toDo.end(), molSims.begin(),
                   calcMolMolSimilarity);
  }
#else
  std::transform(toDo.begin(), toDo.end(), molSims.begin(),
                 calcMolMolSimilarity);
#endif
  for (const auto &cn : molSims) {
    proxGraph[cn.d_mol1Num][cn.d_mol2Num] =
        proxGraph[cn.d_mol2Num][cn.d_mol1Num] = cn;
  }
  return proxGraph;
}

// Split the proximity graph into its disconnected components,
// returning vectors of the molecule numbers of the disconnected
// graphs.
std::vector<std::vector<unsigned int>> disconnectProximityGraphs(
    std::vector<std::vector<ClusNode>> &proxGraph) {
  std::vector<std::vector<unsigned int>> subGraphs;
  std::vector<bool> done(proxGraph.size(), false);
  auto nextStart = std::find(done.begin(), done.end(), false);
  while (nextStart != done.end()) {
    std::list<unsigned int> nodes;
    std::list<unsigned int> toDo(1, std::distance(done.begin(), nextStart));
    while (!toDo.empty()) {
      auto nextNode = toDo.front();
      toDo.pop_front();
      if (!done[nextNode]) {
        nodes.push_back(nextNode);
      }
      done[nextNode] = true;
      for (size_t i = 0; i < proxGraph.size(); ++i) {
        if (!done[i] && proxGraph[nextNode][i].d_res) {
          toDo.push_back(i);
          nodes.push_back(i);
          done[i] = true;
        }
      }
    }
    nodes.sort();
    subGraphs.push_back(std::vector(nodes.begin(), nodes.end()));
    nextStart = std::find(done.begin(), done.end(), false);
  }
  return subGraphs;
}

// Calculate G_{ij} for the molecule.  p is the number of bonds that
// a fragment must exceed for it to be counted in the formula.
double g_ij(const std::shared_ptr<ROMol> &mol, double a, double b,
            unsigned int p) {
  auto molFrags = MolOps::getMolFrags(*mol, false);
  int numBigFrags = 0;
  for (const auto &mf : molFrags) {
    if (mf->getNumBonds() > p) {
      ++numBigFrags;
    }
  }
  numBigFrags = numBigFrags == 0 ? molFrags.size() : numBigFrags;
  double g = mol->getNumAtoms();
  g += b * (1.0 - a * (numBigFrags - 1)) * mol->getNumBonds();
  return g;
}

std::vector<std::vector<unsigned int>> makeSubClusters(
    const std::vector<ClusNode> &nbors, const RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> subClusters;

  std::vector<const ClusNode *> tmpNbors;
  for (const auto &n : nbors) {
    tmpNbors.push_back(&n);
  }

  while (!tmpNbors.empty()) {
    subClusters.push_back(std::vector<unsigned int>{
        tmpNbors.front()->d_mol1Num, tmpNbors.front()->d_mol2Num});
    auto m1 = tmpNbors.front()->d_res->getMcesMol();
    auto g_12 = g_ij(m1, clusOpts.a, clusOpts.b, clusOpts.minFragSize);
    for (size_t i = 1; i < tmpNbors.size(); ++i) {
      auto m2 = tmpNbors[i]->d_res->getMcesMol();
      auto g_13 = g_ij(m2, clusOpts.a, clusOpts.b, clusOpts.minFragSize);

      auto results = RDKit::RascalMCES::rascalMCES(*m1, *m2);
      if (results.empty() || results.front().getBondMatches().empty()) {
        continue;
      }
      auto res = results.front();
      auto g_12_13 =
          g_ij(res.getMcesMol(), clusOpts.a, clusOpts.b, clusOpts.minFragSize);
      double sim = g_12_13 / std::min(g_12, g_13);
      if (sim > clusOpts.minIntraClusterSim) {
        subClusters.back().push_back(tmpNbors[i]->d_mol2Num);
        subClusters.back().push_back(tmpNbors[i]->d_mol1Num);
        tmpNbors[i] = nullptr;
      }
    }
    tmpNbors.front() = nullptr;
    tmpNbors.erase(std::remove(tmpNbors.begin(), tmpNbors.end(), nullptr),
                   tmpNbors.end());
    std::sort(subClusters.back().begin(), subClusters.back().end());
    subClusters.back().erase(
        std::unique(subClusters.back().begin(), subClusters.back().end()),
        subClusters.back().end());
  }
  return subClusters;
}

std::vector<std::vector<unsigned int>> formInitialClusters(
    const std::vector<unsigned int> &subGraph,
    const std::vector<std::vector<ClusNode>> &proxGraph,
    const RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> clusters;
  if (subGraph.size() < 2) {
    return clusters;
  }
  for (auto i : subGraph) {
    std::vector<ClusNode> nbors;
    for (auto j : subGraph) {
      if (proxGraph[i][j].d_res) {
        nbors.push_back(proxGraph[i][j]);
      }
    }
    std::sort(nbors.begin(), nbors.end(),
              [](const ClusNode &c1, const ClusNode &c2) -> bool {
                return c1.d_sim > c2.d_sim;
              });
    if (!nbors.empty()) {
      auto subClusters = makeSubClusters(nbors, clusOpts);
      clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());
    }
  }
  std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<unsigned int> &c1,
               const std::vector<unsigned int> &c2) -> bool {
              if (c1.size() == c2.size()) {
                return c1.front() < c2.front();
              } else {
                return c1.size() > c2.size();
              }
            });
  clusters.erase(std::unique(clusters.begin(), clusters.end()), clusters.end());
  return clusters;
}

std::vector<std::vector<unsigned int>> mergeClusters(
    const std::vector<std::vector<unsigned int>> &clusters,
    const RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> outClusters(clusters);

  if (outClusters.size() < 2) {
    return outClusters;
  }

  for (size_t i = 0; i < outClusters.size() - 1; ++i) {
    for (size_t j = i + 1; j < outClusters.size(); ++j) {
      std::vector<int> inCommon;
      std::set_intersection(outClusters[i].begin(), outClusters[i].end(),
                            outClusters[j].begin(), outClusters[j].end(),
                            std::back_inserter(inCommon));
      double s =
          double(inCommon.size()) / std::min(double(outClusters[i].size()),
                                             double(outClusters[j].size()));
      if (s > clusOpts.clusterMergeSim) {
        outClusters[i].insert(outClusters[i].end(), outClusters[j].begin(),
                              outClusters[j].end());
        outClusters[j].clear();
        std::sort(outClusters[i].begin(), outClusters[i].end());
        outClusters[i].erase(
            std::unique(outClusters[i].begin(), outClusters[i].end()),
            outClusters[i].end());
      }
    }
    outClusters.erase(
        std::remove_if(outClusters.begin(), outClusters.end(),
                       [](const std::vector<unsigned int> &c) -> bool {
                         return c.empty();
                       }),
        outClusters.end());
  }

  return outClusters;
}

void sortClusterMembersByMeanSim(
    const std::vector<std::vector<ClusNode>> &proxGraph,
    std::vector<std::vector<unsigned int>> &clusters) {
  for (auto &clus : clusters) {
    std::vector<std::pair<unsigned int, double>> clusSims;
    for (unsigned int i = 0U; i < clus.size(); ++i) {
      double totSim = 0.0;
      for (unsigned int j = 0U; j < clus.size(); ++j) {
        if (i != j) {
          totSim += proxGraph[clus[i]][clus[j]].d_sim;
        }
      }
      clusSims.push_back({clus[i], totSim / (clus.size() - 1)});
    }
    std::sort(clusSims.begin(), clusSims.end(),
              [](const std::pair<unsigned int, double> &p1,
                 const std::pair<unsigned int, double> &p2) -> bool {
                return p1.second > p2.second;
              });
    std::transform(
        clusSims.begin(), clusSims.end(), clus.begin(),
        [](const std::pair<unsigned int, double> &p) -> unsigned int {
          return p.first;
        });
  }
}

std::vector<std::vector<unsigned int>> makeClusters(
    const std::vector<std::vector<unsigned int>> &subGraphs,
    const std::vector<std::vector<ClusNode>> &proxGraph,
    const RascalClusterOptions &clusOpts) {
  std::vector<std::vector<unsigned int>> clusters;
  for (const auto &sg : subGraphs) {
    auto theseClusters = formInitialClusters(sg, proxGraph, clusOpts);
    auto mergedClusters = mergeClusters(theseClusters, clusOpts);
    clusters.insert(clusters.end(), mergedClusters.begin(),
                    mergedClusters.end());
  }
  std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<unsigned int> &c1,
               const std::vector<unsigned int> &c2) -> bool {
              return c1.size() > c2.size();
            });
  return clusters;
}

std::vector<unsigned int> collectSingletons(
    const std::vector<std::vector<ClusNode>> &proxGraph) {
  std::vector<unsigned int> singletons;
  for (size_t i = 0; i < proxGraph.size(); ++i) {
    bool single = true;
    for (const auto &cn : proxGraph[i]) {
      if (cn.d_res) {
        single = false;
        break;
      }
    }
    if (single) {
      singletons.push_back(i);
    }
  }
  return singletons;
}
}  // namespace details

std::vector<std::vector<unsigned int>> rascalCluster(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts) {
  auto proxGraph = details::buildProximityGraph(mols, clusOpts);
  auto subGraphs = details::disconnectProximityGraphs(proxGraph);
  auto clusters = details::makeClusters(subGraphs, proxGraph, clusOpts);
  auto singletons = details::collectSingletons(proxGraph);
  clusters.push_back(singletons);
  details::sortClusterMembersByMeanSim(proxGraph, clusters);
  return clusters;
}

}  // namespace RascalMCES
}  // namespace RDKit
