
//
//  2019, Daniel Probst, Reymond Group @ University of Bern
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <stdint.h>
#include <math.h> 
#include <set>

#include <RDGeneral/types.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "MHFP.h"

namespace RDKit {
namespace MHFPFingerprints {

MHFPEncoder::MHFPEncoder(unsigned int n_permutations,
             unsigned int seed)
  : n_permutations_(n_permutations)
  , seed_(seed)
  , perms_a_(n_permutations, 0)
  , perms_b_(n_permutations, 0)
{
  std::mt19937 rand;
  rand.seed(seed_);

  std::uniform_int_distribution<std::mt19937::result_type> dist_a(1, max_hash_);
  std::uniform_int_distribution<std::mt19937::result_type> dist_b(0, max_hash_);

  for (unsigned int i = 0; i < n_permutations_; i++) {
    uint32_t a = dist_a(rand);
    uint32_t b = dist_b(rand);

    while (std::find(perms_a_.begin(), perms_a_.end(), a) != perms_a_.end()) {
      a = dist_a(rand);
    }

    while (std::find(perms_b_.begin(), perms_b_.end(), b) != perms_b_.end()) {
      b = dist_a(rand);
    }

    perms_a_[i] = a;
    perms_b_[i] = b;
  }
}

std::vector<uint32_t>
MHFPEncoder::FromStringArray(const std::vector<std::string>& vec)
{
  std::vector<uint32_t> mh(n_permutations_, max_hash_);

  for (uint32_t i = 0; i < vec.size(); i++) {
    auto hashi = FNV::hash(vec[i]);
    for (size_t j = 0; j < n_permutations_; j++){
      uint32_t tmp = (FastMod((perms_a_[j] * hashi + perms_b_[j]), prime_)) & max_hash_;
      mh[j] = std::min(tmp, mh[j]);
    }
  }

  return mh;
}

std::vector<uint32_t>
MHFPEncoder::FromArray(const std::vector<uint32_t>& vec)
{
  std::vector<uint32_t> mh(n_permutations_, max_hash_);

  for (uint32_t i = 0; i < vec.size(); i++) {
    for (size_t j = 0; j < n_permutations_; j++){
      uint32_t tmp = (FastMod((perms_a_[j] * vec[i] + perms_b_[j]), prime_)) & max_hash_;
      mh[j] = std::min(tmp, mh[j]);
    }
  }

  return mh;
}

std::vector<std::string>
MHFPEncoder::CreateShingling(const ROMol& mol, 
               unsigned char radius,
               bool rings,
               bool isomeric,
               bool kekulize,
               unsigned char min_radius) {
  std::vector<std::string> shingling;

  if (rings) {
    const VECT_INT_VECT bonds_vect = mol.getRingInfo()->bondRings();

    for (size_t i = 0; i < bonds_vect.size(); i++) {
      std::unique_ptr<ROMol> m(Subgraphs::pathToSubmol(mol, bonds_vect[i]));
      shingling.emplace_back(MolToSmiles(*m));
    }
  }

  unsigned char min_radius_internal = min_radius;

  if (!min_radius) {
    for (auto atom : mol.atoms()) {
      bool do_kekule = false;
      const RDKit::Bond* bond_in = nullptr;
      bool all_hs_explicit = false;
      bool isomeric_smiles = true;

      shingling.emplace_back(SmilesWrite::GetAtomSmiles(atom, do_kekule, bond_in, all_hs_explicit, isomeric_smiles));
    }

    min_radius_internal++;
  }

  uint32_t index = 0;
  for (auto atom : mol.atoms()) {
    for (unsigned char r = min_radius_internal; r < radius + 1; r++) {
      const PATH_TYPE path = findAtomEnvironmentOfRadiusN(mol, r, index);
      INT_MAP_INT amap;
      bool use_query = false;
      std::unique_ptr<ROMol> submol(Subgraphs::pathToSubmol(mol, path, use_query, amap));

      if (amap.find(index) == amap.end())
        continue;

      std::string smiles = MolToSmiles(*submol, isomeric, kekulize, amap[index]);
      
      if (!smiles.empty())
        shingling.emplace_back(smiles);
    }

    index++;
  }

  return shingling;
}

std::vector<std::string>
MHFPEncoder::CreateShingling(const std::string& smiles, 
               unsigned char radius,
               bool rings,
               bool isomeric,
               bool kekulize,
               unsigned char min_radius) {
  std::unique_ptr<ROMol> m(SmilesToMol(smiles));
  PRECONDITION(m, "could not parse smiles");
  return CreateShingling(*m, radius, rings, isomeric, kekulize, min_radius);
}

std::vector<uint32_t>
MHFPEncoder::Encode(ROMol& mol, 
          unsigned char radius,
          bool rings,
          bool isomeric,
          bool kekulize,
          unsigned char min_radius) {
  return FromStringArray(CreateShingling(mol, radius, rings, isomeric, kekulize, min_radius));                        
}

std::vector<std::vector<uint32_t>>
MHFPEncoder::Encode(std::vector<ROMol>& mols, 
          unsigned char radius,
          bool rings,
          bool isomeric,
          bool kekulize,
          unsigned char min_radius) {
  size_t n = mols.size();
  std::vector<std::vector<uint32_t>> results(n);
  
  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    results[i] = FromStringArray(CreateShingling(mols[i], radius, rings, isomeric, kekulize, min_radius));
  }

  return results;
}

std::vector<uint32_t>
MHFPEncoder::Encode(std::string& smiles, 
          unsigned char radius,
          bool rings,
          bool isomeric,
          bool kekulize,
          unsigned char min_radius) {
  return FromStringArray(CreateShingling(smiles, radius, rings, isomeric, kekulize, min_radius));                        
}

// Someone has to come up with a plural for smiles... smiless, smileses?
std::vector<std::vector<uint32_t>>
MHFPEncoder::Encode(std::vector<std::string>& smileses, 
          unsigned char radius,
          bool rings,
          bool isomeric,
          bool kekulize,
          unsigned char min_radius) {
  size_t n = smileses.size();
  std::vector<std::vector<uint32_t>> results(n);
  
  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    results[i] = FromStringArray(CreateShingling(smileses[i], radius, rings, isomeric, kekulize, min_radius));
  }

  return results;
}

ExplicitBitVect
MHFPEncoder::EncodeSECFP(ROMol& mol, 
             unsigned char radius,
             bool rings,
             bool isomeric,
             bool kekulize,
             unsigned char min_radius,
             size_t length) {
  return Fold(HashShingling(CreateShingling(mol, radius, rings, isomeric, kekulize, min_radius)), length);                        
}

std::vector<ExplicitBitVect>
MHFPEncoder::EncodeSECFP(std::vector<ROMol>& mols, 
             unsigned char radius,
             bool rings,
             bool isomeric,
             bool kekulize,
             unsigned char min_radius,
             size_t length) {
  size_t n = mols.size();
  std::vector<ExplicitBitVect> results(n);
  
  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    results[i] = Fold(HashShingling(CreateShingling(mols[i], radius, rings, isomeric, kekulize, min_radius)), length);
  }

  return results;
}

ExplicitBitVect
MHFPEncoder::EncodeSECFP(std::string& smiles, 
             unsigned char radius,
             bool rings,
             bool isomeric,
             bool kekulize,
             unsigned char min_radius,
             size_t length) {
  return Fold(HashShingling(CreateShingling(smiles, radius, rings, isomeric, kekulize, min_radius)), length);                        
}

std::vector<ExplicitBitVect>
MHFPEncoder::EncodeSECFP(std::vector<std::string>& smileses, 
             unsigned char radius,
             bool rings,
             bool isomeric,
             bool kekulize,
             unsigned char min_radius,
             size_t length) {
  size_t n = smileses.size();
  std::vector<ExplicitBitVect> results(n);
  
  #pragma omp parallel for
  for (size_t i = 0; i < n; i++) {
    results[i] = Fold(HashShingling(CreateShingling(smileses[i], radius, rings, isomeric, kekulize, min_radius)), length);
  }

  return results;
}

}  // namespace MHFPFingerprints
}  // namespace RDKit