//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SCAFFOLDNETWORK_H
#define RD_SCAFFOLDNETWORK_H

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <memory>
#include <iostream>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/version.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {
class ROMol;
class ChemicalReaction;

namespace ScaffoldNetwork {

struct RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetworkParams {
  bool includeGenericScaffolds =
      true;  ///< include scaffolds with all atoms replaced by dummies
  bool includeGenericBondScaffolds =
      false;  ///< include scaffolds with all bonds replaced by single bonds
  bool includeScaffoldsWithoutAttachments =
      true;  ///< remove attachment points from scaffolds and include the result
  bool includeScaffoldsWithAttachments =
      true;  ///< Include the version of the scaffold with attachment points
  bool keepOnlyFirstFragment =
      true;  ///<  keep only the first fragment from the bond breaking rule
  bool pruneBeforeFragmenting =
      true;  ///<  Do a pruning/flattening step before starting fragmenting
  bool flattenIsotopes = true;  ///< remove isotopes when flattening
  bool flattenChirality =
      true;  ///< remove chirality and bond stereo when flattening
  bool flattenKeepLargest =
      true;  ///< keep only the largest fragment when doing flattening
  bool collectMolCounts = true;  ///< keep track of the number of molecules each
                                 ///< scaffold was reached from

  std::vector<std::shared_ptr<ChemicalReaction>>
      bondBreakersRxns;  ///< the reaction(s) used to fragment. Should expect a
                         ///< single reactant and produce two products
  ScaffoldNetworkParams()
      : ScaffoldNetworkParams{
            {"[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]"}} {};
  ScaffoldNetworkParams(const std::vector<std::string> &bondBreakersSmarts);
};

enum class RDKIT_SCAFFOLDNETWORK_EXPORT EdgeType {
  Fragment = 1,     ///< molecule -> fragment
  Generic = 2,      ///< molecule -> generic molecule (all atoms are dummies)
  GenericBond = 3,  ///< molecule -> generic bond molecule (all bonds single)
  RemoveAttachment = 4,  ///< molecule -> molecule with no attachment points
  Initialize = 5         ///< molecule -> flattened molecule
};

struct RDKIT_SCAFFOLDNETWORK_EXPORT NetworkEdge {
  size_t beginIdx;
  size_t endIdx;
  EdgeType type;
  NetworkEdge() : beginIdx(0), endIdx(0), type(EdgeType::Initialize){};
  NetworkEdge(size_t bi, size_t ei, EdgeType typ)
      : beginIdx(bi), endIdx(ei), type(typ){};
  bool operator==(const RDKit::ScaffoldNetwork::NetworkEdge &o) const {
    return (beginIdx == o.beginIdx) && (endIdx == o.endIdx) && (type == o.type);
  }
  bool operator!=(const RDKit::ScaffoldNetwork::NetworkEdge &o) const {
    return (beginIdx != o.beginIdx) || (endIdx != o.endIdx) || (type != o.type);
  }
#ifdef RDK_USE_BOOST_SERIALIZATION
 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &beginIdx;
    ar &endIdx;
    ar &type;
  }
#endif
};

struct RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetwork {
  std::vector<std::string> nodes;  ///< SMILES for the scaffolds
  std::vector<unsigned>
      counts;  ///< number of times each scaffold was encountered
  std::vector<unsigned>
      molCounts;  ///< number of molecules each scaffold was found in
  std::vector<NetworkEdge> edges;  ///< edges in the network
  ScaffoldNetwork(){};
#ifdef RDK_USE_BOOST_SERIALIZATION
  ScaffoldNetwork(const std::string &pkl) {
    std::stringstream iss(pkl);
    boost::archive::text_iarchive ia(iss);
    ia >> *this;
  }

 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);
    ar &nodes;
    ar &counts;
    if (version > 0) {
      ar &molCounts;
    }
    ar &edges;
  }
#endif
};

//! update an existing ScaffoldNetwork using a set of molecules
template <typename T>
void updateScaffoldNetwork(const T &mols, ScaffoldNetwork &network,
                           const ScaffoldNetworkParams &params);

//! create a new ScaffoldNetwork for a set of molecules
template <typename T>
ScaffoldNetwork createScaffoldNetwork(const T &mols,
                                      const ScaffoldNetworkParams &params) {
  ScaffoldNetwork res;
  updateScaffoldNetwork(mols, res, params);
  return res;
}
//! allows nodes to output nicely as strings
inline std::ostream &operator<<(std::ostream &ostr,
                                const RDKit::ScaffoldNetwork::EdgeType &e) {
  switch (e) {
    case RDKit::ScaffoldNetwork::EdgeType::Fragment:
      ostr << "Fragment";
      break;
    case RDKit::ScaffoldNetwork::EdgeType::Generic:
      ostr << "Generic";
      break;
    case RDKit::ScaffoldNetwork::EdgeType::GenericBond:
      ostr << "GenericBond";
      break;
    case RDKit::ScaffoldNetwork::EdgeType::RemoveAttachment:
      ostr << "RemoveAttachment";
      break;
    case RDKit::ScaffoldNetwork::EdgeType::Initialize:
      ostr << "Initialize";
      break;
    default:
      ostr << "UNKNOWN";
      break;
  }
  return ostr;
}
//! allows edges to output nicely as strings
inline std::ostream &operator<<(std::ostream &ostr,
                                const RDKit::ScaffoldNetwork::NetworkEdge &e) {
  ostr << "NetworkEdge( " << e.beginIdx << "->" << e.endIdx
       << ", type:" << e.type << " )";
  return ostr;
}

//! returns parameters for constructing scaffold networks using BRICS
//! fragmentation
RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetworkParams getBRICSNetworkParams();

}  // namespace ScaffoldNetwork
}  // namespace RDKit

namespace boost {
namespace serialization {
template <>
struct version<RDKit::ScaffoldNetwork::ScaffoldNetwork> {
  BOOST_STATIC_CONSTANT(int, value = 1);
};
}  // namespace serialization
}  // namespace boost

#endif
