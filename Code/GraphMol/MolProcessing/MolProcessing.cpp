//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MolProcessing.h"

namespace RDKit {
namespace MolProcessing {

namespace details {
GeneralMolSupplier::SupplierOptions defaultSupplierOptions;
inline std::unique_ptr<FileParsers::MolSupplier> getSupplier(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static std::once_flag flag;
  std::call_once(flag, []() { defaultSupplierOptions.numWriterThreads = 4; });
#endif
  return GeneralMolSupplier::getSupplier(fileName, options);
}
}  // namespace details

namespace {
#ifdef RDK_BUILD_THREADSAFE_SSS
inline std::mutex &get_fp_mutex() {
  // create on demand
  static std::mutex _mutex;
  return _mutex;
}

template <typename T>
std::vector<std::unique_ptr<T>> mtWorker(
    v2::FileParsers::MultithreadedMolSupplier *suppl,
    std::function<T *(RWMol &)> func) {
  PRECONDITION(suppl, "no supplier");
  std::map<unsigned int, std::unique_ptr<T>> accum;

  auto workerfunc = [&](RWMol &mol, const std::string &,
                        unsigned int recordId) {
    auto item = func(mol);
    {
      std::lock_guard<std::mutex> lock(get_fp_mutex());
      accum[recordId].reset(item);
    }
  };
  suppl->setWriteCallback(workerfunc);
  // loop over the supplier to make sure we read everything
  while (!suppl->atEnd()) {
    auto mol = suppl->next();
  }
  // convert the map to a vector and get the results in the input order
  auto maxv = 0u;
  for (const auto &pr : accum) {
    maxv = std::max(maxv, pr.first);
  }
  std::vector<std::unique_ptr<T>> results(maxv);
  for (auto &pr : accum) {
    results[pr.first - 1] = std::move(pr.second);
  }
  return results;
}
#endif

template <typename T>
std::vector<std::unique_ptr<T>> worker(v2::FileParsers::MolSupplier *suppl,
                                       std::function<T *(RWMol &)> func) {
  PRECONDITION(suppl, "no supplier");
  // if we are using a multi-threaded supplier then we can register a write
  // callback to do our processing multi-threaded too
#ifdef RDK_BUILD_THREADSAFE_SSS
  auto tsuppl =
      dynamic_cast<v2::FileParsers::MultithreadedMolSupplier *>(suppl);
  if (tsuppl) {
    return mtWorker(tsuppl, func);
  } else {
#else
  {
#endif
    // otherwise we just loop through the molecules
    std::vector<std::unique_ptr<T>> results;
    while (!suppl->atEnd()) {
      auto mol = suppl->next();
      if (mol) {
        auto fp = func(*mol);
        results.emplace_back(fp);
      } else if (!suppl->atEnd()) {
        results.emplace_back(nullptr);
      }
    }
    return results;
  }
}
}  // namespace

//! \brief Get fingerprints for all of the molecules in a file
/*!
   \param fileName the name of the file to read
   \param options options controlling how the file is read, if not provided
           four threads will be used whegn reading the file
   \param generator the fingerprint generator to use, if not provided,
           Morgan fingerprints with radius of 3 will be used.

   \return an ExplicitBitVect,bitset pair, the first containing the
           fingerprints and the second a bitset indicating which molecules were
           successfully read
*/
template <typename OutputType>
std::vector<std::unique_ptr<ExplicitBitVect>> getFingerprintsForMolsInFile(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options,
    FingerprintGenerator<OutputType> *generator) {
  auto suppl = details::getSupplier(fileName, options);

  std::unique_ptr<FingerprintGenerator<OutputType>> morgan;
  if (generator == nullptr) {
    morgan.reset(MorganFingerprint::getMorganGenerator<OutputType>(3));
    generator = morgan.get();
  }
  std::function<ExplicitBitVect *(RWMol &)> func = [&](RWMol &mol) {
    return generator->getFingerprint(mol);
  };
  auto results = worker(suppl.get(), func);
  return results;
}

template RDKIT_MOLPROCESSING_EXPORT
    std::vector<std::unique_ptr<ExplicitBitVect>>
    getFingerprintsForMolsInFile(
        const std::string &fileName,
        const GeneralMolSupplier::SupplierOptions &options,
        FingerprintGenerator<std::uint32_t> *generator);
template RDKIT_MOLPROCESSING_EXPORT
    std::vector<std::unique_ptr<ExplicitBitVect>>
    getFingerprintsForMolsInFile(
        const std::string &fileName,
        const GeneralMolSupplier::SupplierOptions &options,
        FingerprintGenerator<std::uint64_t> *generator);

}  // namespace MolProcessing
}  // namespace RDKit
