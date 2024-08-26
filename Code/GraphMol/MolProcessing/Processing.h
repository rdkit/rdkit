//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLPROCCESING_H
#define RD_MOLPROCCESING_H

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/GeneralFileReader.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

namespace RDKit {
namespace MolProccesing {

auto defaultSupplierOptions = GeneralMolSupplier::SupplierOptions();

namespace details {
inline std::unique_ptr<FileParsers::MolSupplier> getSupplier(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options) {
  static bool firstCall = true;
  if (firstCall) {
    defaultSupplierOptions.numWriterThreads = 4;
    firstCall = false;
  }
  return GeneralMolSupplier::getSupplier(fileName, options);
}
}  // namespace details

//! \brief Get fingerprints for all of the molecules in a file
/*!
   \param fileName the name of the file to read
   \param options options controlling how the file is read, if not provided
           four threads will be used when reading the file
   \param generator the fingerprint generator to use, if not provided,
           Morgan fingerprints with radius of 3 will be used.

   \return an ExplicitBitVect,bitset pair, the first containing the
           fingerprints and the second a bitset indicating which molecules were
           successfully read
*/
template <typename OutputType = std::uint32_t>
std::pair<std::vector<std::unique_ptr<ExplicitBitVect>>,
          boost::dynamic_bitset<>>
getFingerprintsForMolsInFile(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options = defaultSupplierOptions,
    FingerprintGenerator<OutputType> *generator = nullptr) {
  auto suppl = details::getSupplier(fileName, options);

  std::unique_ptr<FingerprintGenerator<OutputType>> morgan;
  if (generator == nullptr) {
    morgan.reset(MorganFingerprint::getMorganGenerator<OutputType>(3));
    generator = morgan.get();
  }
  std::vector<std::unique_ptr<ExplicitBitVect>> fingerprints;
  boost::dynamic_bitset<> passed;
  while (!suppl->atEnd()) {
    auto mol = suppl->next();
    if (mol) {
      fingerprints.emplace_back(generator->getFingerprint(*mol));
      passed.push_back(true);
    } else {
      passed.push_back(false);
    }
  }
  return std::make_pair(std::move(fingerprints), passed);
}

}  // namespace MolProccesing
}  // namespace RDKit
#endif
