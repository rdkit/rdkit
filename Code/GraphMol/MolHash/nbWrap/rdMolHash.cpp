//
//  Copyright (C) 2020-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolHash/MolHash.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

NB_MODULE(rdMolHash, m) {
  m.doc() = "Module containing functions to generate hashes for molecules";

  nb::enum_<MolHash::HashFunction>(m, "HashFunction")
      .value("AnonymousGraph", MolHash::HashFunction::AnonymousGraph)
      .value("ElementGraph", MolHash::HashFunction::ElementGraph)
      .value("CanonicalSmiles", MolHash::HashFunction::CanonicalSmiles)
      .value("MurckoScaffold", MolHash::HashFunction::MurckoScaffold)
      .value("ExtendedMurcko", MolHash::HashFunction::ExtendedMurcko)
      .value("MolFormula", MolHash::HashFunction::MolFormula)
      .value("AtomBondCounts", MolHash::HashFunction::AtomBondCounts)
      .value("DegreeVector", MolHash::HashFunction::DegreeVector)
      .value("Mesomer", MolHash::HashFunction::Mesomer)
      .value("HetAtomTautomer", MolHash::HashFunction::HetAtomTautomer)
      .value("HetAtomProtomer", MolHash::HashFunction::HetAtomProtomer)
      .value("RedoxPair", MolHash::HashFunction::RedoxPair)
      .value("Regioisomer", MolHash::HashFunction::Regioisomer)
      .value("NetCharge", MolHash::HashFunction::NetCharge)
      .value("SmallWorldIndexBR", MolHash::HashFunction::SmallWorldIndexBR)
      .value("SmallWorldIndexBRL", MolHash::HashFunction::SmallWorldIndexBRL)
      .value("ArthorSubstructureOrder",
             MolHash::HashFunction::ArthorSubstructureOrder)
      .value("HetAtomTautomerv2", MolHash::HashFunction::HetAtomTautomerv2)
      .value("HetAtomProtomerv2", MolHash::HashFunction::HetAtomProtomerv2)
      .export_values();

  m.def(
      "MolHash",
      [](const ROMol &mol, MolHash::HashFunction func, bool useCXSmiles,
         unsigned cxFlagsToSkip) {
        RWMol cpy(mol);
        return MolHash::MolHash(&cpy, func, useCXSmiles, cxFlagsToSkip);
      },
      "mol"_a, "func"_a, "useCxSmiles"_a = false, "cxFlagsToSkip"_a = 0u,
      R"DOC(Generate a hash for a molecule. The func argument determines which hash is generated.)DOC");
}
