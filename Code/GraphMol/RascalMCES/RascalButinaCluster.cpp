//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains an implementation of Butina clustering
// (Butina JCICS 39 747-750 (1999)) using the RascalMCES
// Johnson similarity metric.  It is largely a transliteration
// of $RDBASE/rdkit/ML/Cluster/Butina.py.

#include <algorithm>
#include <iterator>
#include <vector>
#include <set>

#include <GraphMol/ROMol.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalDetails.h>

namespace RDKit {

namespace RascalMCES {
namespace details {
std::vector<std::vector<unsigned int>> buildNborLists(
    const std::vector<std::vector<ClusNode>> &proxGraph) {
  std::vector<std::vector<unsigned int>> nborLists;
  for (size_t i = 0; i < proxGraph.size(); ++i) {
    std::vector<std::pair<unsigned int, double>> tmpList;
    for (const auto &cn : proxGraph[i]) {
      if (cn.d_res) {
        if (i == cn.d_mol1Num) {
          tmpList.push_back({cn.d_mol2Num, cn.d_sim});
        } else {
          tmpList.push_back({cn.d_mol1Num, cn.d_sim});
        }
      }
    }
    std::sort(tmpList.begin(), tmpList.end(),
              [](const std::pair<unsigned int, double> &p1,
                 const std::pair<unsigned int, double> &p2) -> bool {
                return p1.second > p2.second;
              });
    std::vector<unsigned int> nborList(tmpList.size() + 1, 0);
    nborList[0] = i;
    std::transform(
        tmpList.begin(), tmpList.end(), nborList.begin() + 1,
        [](const std::pair<unsigned int, double> &p) -> unsigned int {
          return p.first;
        });
    nborLists.push_back(nborList);
  }
  std::sort(nborLists.begin(), nborLists.end(),
            [](const std::vector<unsigned int> &nl1,
               const std::vector<unsigned int> &nl2) -> bool {
              if (nl1.size() == nl2.size()) {
                return nl1 > nl2;
              } else {
                return nl1.size() > nl2.size();
              }
            });
  return nborLists;
}

// This function destroys nborLists.
std::vector<std::vector<unsigned int>> formClusters(
    std::vector<std::vector<unsigned int>> &nborLists) {
  std::vector<std::vector<unsigned int>> clusters;

  while (!nborLists.empty()) {
    clusters.push_back(nborLists.front());
    std::set<unsigned int> inNborList(nborLists.front().begin(),
                                      nborLists.front().end());
    nborLists.front().clear();
    for (auto &nborList : nborLists) {
      for (auto &n : nborList) {
        if (inNborList.find(n) != inNborList.end()) {
          n = std::numeric_limits<unsigned int>::max();
        }
      }
      nborList.erase(std::remove(nborList.begin(), nborList.end(),
                                 std::numeric_limits<unsigned int>::max()),
                     nborList.end());
    }
    nborLists.erase(
        std::remove_if(nborLists.begin(), nborLists.end(),
                       [](const std::vector<unsigned int> &nl) -> bool {
                         return nl.empty();
                       }),
        nborLists.end());
    std::sort(nborLists.begin(), nborLists.end(),
              [](const std::vector<unsigned int> &nl1,
                 const std::vector<unsigned int> &nl2) -> bool {
                if (nl1.size() == nl2.size()) {
                  return nl1 > nl2;
                } else {
                  return nl1.size() > nl2.size();
                }
              });
  }
  return clusters;
}

}  // namespace details
std::vector<std::vector<unsigned int>> rascalButinaCluster(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts) {
  auto proxGraph = details::buildProximityGraph(mols, clusOpts);
  auto nborLists = details::buildNborLists(proxGraph);
  auto clusters = details::formClusters(nborLists);
  return clusters;
}
}  // namespace RascalMCES
}  // namespace RDKit
