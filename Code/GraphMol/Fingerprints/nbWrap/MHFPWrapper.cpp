//
//  2019, Daniel Probst, Reymond Group @ University of Bern
//  Copyright (C) 2019-2026 and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <GraphMol/Fingerprints/MHFP.h>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;
using RDKit::MHFPFingerprints::MHFPEncoder;

namespace RDKit {
namespace MHFPWrapper {

template <typename T>
std::vector<T> ListToVector(nb::object obj) {
  std::vector<T> result;
  for (auto item : obj) {
    result.push_back(nb::cast<T>(item));
  }
  return result;
}

NB_MODULE(rdMHFPFingerprint, m) {
  nb::class_<MHFPEncoder>(m, "MHFPEncoder")
      .def(nb::init<unsigned int, unsigned int>(), "n_permutations"_a = 2048,
           "seed"_a = 42)
      .def(
          "FromStringArray",
          [](MHFPEncoder *enc, nb::object vec) {
            auto v = ListToVector<std::string>(vec);
            return enc->FromStringArray(v);
          },
          "vec"_a, "Creates a MHFP vector from a list of arbitrary strings.")
      .def(
          "FromArray",
          [](MHFPEncoder *enc, nb::object vec) {
            auto v = ListToVector<uint32_t>(vec);
            return enc->FromArray(v);
          },
          "vec"_a, "Creates a MHFP vector from a list of unsigned integers.")
      .def(
          "CreateShinglingFromSmiles",
          [](MHFPEncoder *enc, std::string smiles, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) {
            return enc->CreateShingling(smiles, radius, rings, isomeric,
                                        kekulize, min_radius);
          },
          "smiles"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          R"DOC(Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.)DOC")
      .def(
          "CreateShinglingFromMol",
          [](MHFPEncoder *enc, const ROMol &mol, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) {
            ROMol molCopy = mol;
            return enc->CreateShingling(molCopy, radius, rings, isomeric,
                                        kekulize, min_radius);
          },
          "mol"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          R"DOC(Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.)DOC")
      .def(
          "EncodeSmiles",
          [](MHFPEncoder *enc, std::string smiles, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) {
            return enc->Encode(smiles, radius, rings, isomeric, kekulize,
                               min_radius);
          },
          "smiles"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          "Creates a MHFP vector from a SMILES string.")
      .def(
          "EncodeMol",
          [](MHFPEncoder *enc, const ROMol &mol, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) {
            ROMol molCopy = mol;
            return enc->Encode(molCopy, radius, rings, isomeric, kekulize,
                               min_radius);
          },
          "mol"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          "Creates a MHFP vector from an RDKit Mol instance.")
      .def(
          "EncodeSmilesBulk",
          [](MHFPEncoder *enc, nb::object smiles, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) -> nb::tuple {
            auto vec = ListToVector<std::string>(smiles);
            auto resVect =
                enc->Encode(vec, radius, rings, isomeric, kekulize, min_radius);
            nb::list resList;
            for (const auto &item : resVect) {
              resList.append(item);
            }
            return nb::tuple(resList);
          },
          "smiles"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          "Creates a MHFP vector from a list of SMILES strings.")
      .def(
          "EncodeMolsBulk",
          [](MHFPEncoder *enc, nb::object mols, unsigned char radius,
             bool rings, bool isomeric, bool kekulize,
             unsigned char min_radius) -> nb::tuple {
            // There are access problems for vector and std::vector<ROMol>.
            // So let's fall back to an inefficient Python list.
            std::vector<ROMol> vec;
            for (auto item : mols) {
              vec.push_back(*nb::cast<const ROMol *>(item));
            }
            auto resVect =
                enc->Encode(vec, radius, rings, isomeric, kekulize, min_radius);
            nb::list resList;
            for (const auto &item : resVect) {
              resList.append(item);
            }
            return nb::tuple(resList);
          },
          "mols"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1,
          "Creates a MHFP vector from a list of RDKit Mol instances.")
      .def(
          "EncodeSECFPSmiles",
          [](MHFPEncoder *enc, std::string smiles, unsigned char radius,
             bool rings, bool isomeric, bool kekulize, unsigned char min_radius,
             size_t length) {
            return enc->EncodeSECFP(smiles, radius, rings, isomeric, kekulize,
                                    min_radius, length);
          },
          "smiles"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1, "length"_a = 2048,
          "Creates a SECFP binary vector from a SMILES string.")
      .def(
          "EncodeSECFPMol",
          [](MHFPEncoder *enc, const ROMol &mol, unsigned char radius,
             bool rings, bool isomeric, bool kekulize, unsigned char min_radius,
             size_t length) {
            ROMol molCopy = mol;
            return enc->EncodeSECFP(molCopy, radius, rings, isomeric, kekulize,
                                    min_radius, length);
          },
          "mol"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1, "length"_a = 2048,
          "Creates a SECFP binary vector from an RDKit Mol instance.")
      .def(
          "EncodeSECFPSmilesBulk",
          [](MHFPEncoder *enc, nb::object smiles, unsigned char radius,
             bool rings, bool isomeric, bool kekulize, unsigned char min_radius,
             size_t length) -> nb::tuple {
            auto vec = ListToVector<std::string>(smiles);
            auto resVect = enc->EncodeSECFP(vec, radius, rings, isomeric,
                                            kekulize, min_radius, length);
            nb::list resList;
            for (const auto &item : resVect) {
              resList.append(item);
            }
            return nb::tuple(resList);
          },
          "smiles"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1, "length"_a = 2048,
          "Creates a SECFP binary vector from a list of SMILES strings.")
      .def(
          "EncodeSECFPMolsBulk",
          [](MHFPEncoder *enc, nb::object mols, unsigned char radius,
             bool rings, bool isomeric, bool kekulize, unsigned char min_radius,
             size_t length) -> nb::tuple {
            std::vector<ROMol> vec;
            for (auto item : mols) {
              vec.push_back(*nb::cast<const ROMol *>(item));
            }
            auto resVect = enc->EncodeSECFP(vec, radius, rings, isomeric,
                                            kekulize, min_radius, length);
            nb::list resList;
            for (const auto &item : resVect) {
              resList.append(item);
            }
            return nb::tuple(resList);
          },
          "mols"_a, "radius"_a = 3, "rings"_a = true, "isomeric"_a = false,
          "kekulize"_a = true, "min_radius"_a = 1, "length"_a = 2048,
          "Creates a SECFP binary vector from a list of RDKit Mol instances.")
      .def_static("Distance", &MHFPEncoder::Distance, "a"_a, "b"_a);
}

}  // namespace MHFPWrapper
}  // namespace RDKit
