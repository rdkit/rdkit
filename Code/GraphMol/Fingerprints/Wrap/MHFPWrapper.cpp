//
//  2019, Daniel Probst, Reymond Group @ University of Bern
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <GraphMol/Fingerprints/MHFP.h>
#include <vector>

namespace python = boost::python;
using RDKit::MHFPFingerprints::MHFPEncoder;

namespace RDKit {
namespace MHFPWrapper {

typedef std::vector<std::vector<uint32_t>> VectMinHashVect;
typedef std::vector<ExplicitBitVect> VectExplicitBitVect;

template <typename T>
std::vector<T> ListToVector(const python::object& obj) {
  return std::vector<T>(python::stl_input_iterator<T>(obj),
                        python::stl_input_iterator<T>());
}

UINT_VECT
FromStringArray(MHFPEncoder* mhfpEnc, python::list& vec) {
  std::vector<std::string> vec_tmp = ListToVector<std::string>(vec);
  return mhfpEnc->FromStringArray(vec_tmp);
}

UINT_VECT
FromArray(MHFPEncoder* mhfpEnc, python::list& vec) {
  std::vector<uint32_t> vec_tmp = ListToVector<uint32_t>(vec);
  return mhfpEnc->FromArray(vec_tmp);
}

STR_VECT
CreateShinglingFromSmiles(MHFPEncoder* mhfpEnc, std::string smiles,
                          unsigned char radius = 3, bool rings = true,
                          bool isomeric = false, bool kekulize = true,
                          unsigned char min_radius = 1) {
  return mhfpEnc->CreateShingling(smiles, radius, rings, isomeric, kekulize,
                                  min_radius);
}

STR_VECT
CreateShinglingFromMol(MHFPEncoder* mhfpEnc, ROMol mol,
                       unsigned char radius = 3, bool rings = true,
                       bool isomeric = false, bool kekulize = true,
                       unsigned char min_radius = 1) {
  return mhfpEnc->CreateShingling(mol, radius, rings, isomeric, kekulize,
                                  min_radius);
}

UINT_VECT
EncodeSmiles(MHFPEncoder* mhfpEnc, std::string smiles, unsigned char radius = 3,
             bool rings = true, bool isomeric = false, bool kekulize = true,
             unsigned char min_radius = 1) {
  return mhfpEnc->Encode(smiles, radius, rings, isomeric, kekulize, min_radius);
}

UINT_VECT
EncodeMol(MHFPEncoder* mhfpEnc, ROMol mol, unsigned char radius = 3,
          bool rings = true, bool isomeric = false, bool kekulize = true,
          unsigned char min_radius = 1) {
  return mhfpEnc->Encode(mol, radius, rings, isomeric, kekulize, min_radius);
}

VectMinHashVect EncodeSmilesBulk(MHFPEncoder* mhfpEnc, python::list& smiles,
                                 unsigned char radius = 3, bool rings = true,
                                 bool isomeric = false, bool kekulize = true,
                                 unsigned char min_radius = 1) {
  std::vector<std::string> vec = ListToVector<std::string>(smiles);
  return mhfpEnc->Encode(vec, radius, rings, isomeric, kekulize, min_radius);
}

// There are access problems for vector_indexing_suite and std::vector<ROMol>.
// So let's fallback to a unefficient Python list.
VectMinHashVect EncodeMolsBulk(MHFPEncoder* mhfpEnc, python::list& mols,
                               unsigned char radius = 3, bool rings = true,
                               bool isomeric = false, bool kekulize = true,
                               unsigned char min_radius = 1) {
  std::vector<ROMol> vec = ListToVector<ROMol>(mols);
  return mhfpEnc->Encode(vec, radius, rings, isomeric, kekulize, min_radius);
}

// SECFP

ExplicitBitVect EncodeSECFPSmiles(MHFPEncoder* mhfpEnc, std::string smiles,
                                  unsigned char radius = 3, bool rings = true,
                                  bool isomeric = false, bool kekulize = true,
                                  unsigned char min_radius = 1,
                                  size_t length = 2048) {
  return mhfpEnc->EncodeSECFP(smiles, radius, rings, isomeric, kekulize,
                              min_radius, length);
}

ExplicitBitVect EncodeSECFPMol(MHFPEncoder* mhfpEnc, ROMol mol,
                               unsigned char radius = 3, bool rings = true,
                               bool isomeric = false, bool kekulize = true,
                               unsigned char min_radius = 1,
                               size_t length = 2048) {
  return mhfpEnc->EncodeSECFP(mol, radius, rings, isomeric, kekulize,
                              min_radius, length);
}

VectExplicitBitVect EncodeSECFPSmilesBulk(
    MHFPEncoder* mhfpEnc, python::list& smiles, unsigned char radius = 3,
    bool rings = true, bool isomeric = false, bool kekulize = true,
    unsigned char min_radius = 1, size_t length = 2048) {
  std::vector<std::string> vec = ListToVector<std::string>(smiles);
  return mhfpEnc->EncodeSECFP(vec, radius, rings, isomeric, kekulize,
                              min_radius, length);
}

VectExplicitBitVect EncodeSECFPMolsBulk(
    MHFPEncoder* mhfpEnc, python::list& mols, unsigned char radius = 3,
    bool rings = true, bool isomeric = false, bool kekulize = true,
    unsigned char min_radius = 1, size_t length = 2048) {
  std::vector<ROMol> vec = ListToVector<ROMol>(mols);
  return mhfpEnc->EncodeSECFP(vec, radius, rings, isomeric, kekulize,
                              min_radius, length);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(CreateShinglingFromSmilesOverloads,
                                CreateShinglingFromSmiles, 2, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(CreateShinglingFromMolOverloads,
                                CreateShinglingFromMol, 2, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSmilesOverloads, EncodeSmiles, 2, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeMolOverloads, EncodeMol, 2, 7)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSmilesBulkOverloads, EncodeSmilesBulk, 2,
                                7)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeMolsBulkOverloads, EncodeMolsBulk, 2, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSECFPSmilesOverloads, EncodeSECFPSmiles,
                                2, 8)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSECFPMolOverloads, EncodeSECFPMol, 2, 8)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSECFPSmilesBulkOverloads,
                                EncodeSECFPSmilesBulk, 2, 8)
BOOST_PYTHON_FUNCTION_OVERLOADS(EncodeSECFPMolsBulkOverloads,
                                EncodeSECFPMolsBulk, 2, 8)

BOOST_PYTHON_MODULE(rdMHFPFingerprint) {
  python::class_<MHFPEncoder>(
      "MHFPEncoder",
      python::init<python::optional<unsigned int, unsigned int>>())
      .def("FromStringArray", &FromStringArray,
           ((python::arg("vec")),
            "Creates a MHFP vector from a list of arbitrary strings."))
      .def("FromArray", &FromArray,
           ((python::arg("vec")),
            "Creates a MHFP vector from a list of unsigned integers."))
      .def("CreateShinglingFromSmiles", CreateShinglingFromSmiles,
           CreateShinglingFromSmilesOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a shingling (a list of circular n-grams / "
               "substructures) from a SMILES string."))
      .def("CreateShinglingFromMol", CreateShinglingFromMol,
           CreateShinglingFromMolOverloads(
               (python::arg("mol"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a shingling (a list of circular n-grams / "
               "substructures) from a RDKit Mol instance."))
      .def("EncodeSmiles", EncodeSmiles,
           EncodeSmilesOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a MHFP vector from a SMILES string."))
      .def("EncodeMol", EncodeMol,
           EncodeMolOverloads(
               (python::arg("mol"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a MHFP vector from an RDKit Mol instance."))
      .def("EncodeSmilesBulk", EncodeSmilesBulk,
           EncodeSmilesBulkOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a MHFP vector from a list of SMILES strings."))
      .def("EncodeMolsBulk", EncodeMolsBulk,
           EncodeMolsBulkOverloads(
               (python::arg("mols"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1),
               "Creates a MHFP vector from a list of RDKit Mol instances."))
      .def("EncodeSECFPSmiles", EncodeSECFPSmiles,
           EncodeSECFPSmilesOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1,
                python::arg("length") = 2048),
               "Creates a SECFP binary vector from a SMILES string."))
      .def("EncodeSECFPMol", EncodeSECFPMol,
           EncodeSECFPMolOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1,
                python::arg("length") = 2048),
               "Creates a SECFP binary vector from an RDKit Mol instance."))
      .def("EncodeSECFPSmilesBulk", EncodeSECFPSmilesBulk,
           EncodeSECFPSmilesBulkOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1,
                python::arg("length") = 2048),
               "Creates a SECFP binary vector from a list of SMILES strings."))
      .def("EncodeSECFPMolsBulk", EncodeSECFPMolsBulk,
           EncodeSECFPMolsBulkOverloads(
               (python::arg("smiles"), python::arg("radius") = 3,
                python::arg("rings") = true, python::arg("isomeric") = false,
                python::arg("kekulize") = false, python::arg("min_radius") = 1,
                python::arg("length") = 2048),
               "Creates a SECFP binary vector from a list of RDKit Mol "
               "instances."))
      .def("Distance", &MHFPEncoder::Distance)
      .staticmethod("Distance");
}

}  // namespace MHFPWrapper
}  // namespace RDKit
