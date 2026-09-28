//
//  Copyright (C) 2019-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include "EHTTools.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::ndarray<nb::numpy, double, nb::ndim<2>> getMatrixProp(const double *mat,
                                                           unsigned int dim1,
                                                           unsigned int dim2) {
  if (!mat) {
    throw nb::value_error("matrix has not been initialized");
  }
  auto *resData = new double[dim1 * dim2];
  memcpy(resData, static_cast<const void *>(mat),
         dim1 * dim2 * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {dim1, dim2},
                                                     owner);
}

nb::ndarray<nb::numpy, double, nb::ndim<1>> getSymmMatrixProp(
    const double *mat, unsigned int sz) {
  if (!mat) {
    throw nb::value_error("matrix has not been initialized");
  }
  size_t n = (size_t)sz * (sz + 1) / 2;
  auto *resData = new double[n];
  memcpy(resData, static_cast<const void *>(mat), n * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<1>>(resData, {n}, owner);
}

nb::ndarray<nb::numpy, double, nb::ndim<1>> getVectorProp(const double *mat,
                                                           unsigned int sz) {
  if (!mat) {
    throw nb::value_error("vector has not been initialized");
  }
  auto *resData = new double[sz];
  memcpy(resData, static_cast<const void *>(mat), sz * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<1>>(resData, {sz}, owner);
}

}  // end of anonymous namespace

NB_MODULE(rdEHTTools, m) {
  m.doc() =
      R"DOC(Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next.)DOC";

  nb::class_<RDKit::EHTTools::EHTResults>(m, "EHTResults")
      .def_ro("numOrbitals", &RDKit::EHTTools::EHTResults::numOrbitals)
      .def_ro("numElectrons", &RDKit::EHTTools::EHTResults::numElectrons)
      .def_ro("fermiEnergy", &RDKit::EHTTools::EHTResults::fermiEnergy)
      .def_ro("totalEnergy", &RDKit::EHTTools::EHTResults::totalEnergy)
      .def(
          "GetReducedChargeMatrix",
          [](RDKit::EHTTools::EHTResults &self) {
            return getMatrixProp(self.reducedChargeMatrix.get(), self.numAtoms,
                                 self.numOrbitals);
          },
          "returns the reduced charge matrix")
      .def(
          "GetReducedOverlapPopulationMatrix",
          [](RDKit::EHTTools::EHTResults &self) {
            return getSymmMatrixProp(
                self.reducedOverlapPopulationMatrix.get(), self.numAtoms);
          },
          "returns the reduced overlap population matrix")
      .def(
          "GetAtomicCharges",
          [](RDKit::EHTTools::EHTResults &self) {
            return getVectorProp(self.atomicCharges.get(), self.numAtoms);
          },
          "returns the calculated atomic charges")
      .def(
          "GetHamiltonian",
          [](RDKit::EHTTools::EHTResults &self) {
            if (!self.hamiltonianMatrix) {
              throw nb::value_error(
                  "Hamiltonian not available, set "
                  "keepOverlapAndHamiltonianMatrices=True "
                  "to preserve it.");
            }
            return getMatrixProp(self.hamiltonianMatrix.get(),
                                 self.numOrbitals, self.numOrbitals);
          },
          "returns the symmetric Hamiltonian matrix")
      .def(
          "GetOverlapMatrix",
          [](RDKit::EHTTools::EHTResults &self) {
            if (!self.overlapMatrix) {
              throw nb::value_error(
                  "Overlap matrix not available, set "
                  "keepOverlapAndHamiltonianMatrices=True "
                  "to preserve it.");
            }
            return getMatrixProp(self.overlapMatrix.get(), self.numOrbitals,
                                 self.numOrbitals);
          },
          "returns the symmetric overlap matrix")
      .def(
          "GetOrbitalEnergies",
          [](RDKit::EHTTools::EHTResults &self) {
            return getVectorProp(self.orbitalEnergies.get(), self.numOrbitals);
          },
          "returns the energies of the molecular orbitals as a vector");

  m.def(
      "RunMol",
      [](const RDKit::ROMol &mol, int confId,
         bool keepOverlapAndHamiltonianMatrices) {
        auto *eRes = new RDKit::EHTTools::EHTResults();
        bool ok = RDKit::EHTTools::runMol(mol, *eRes, confId,
                                          keepOverlapAndHamiltonianMatrices);
        return nb::make_tuple(
            ok, nb::cast(eRes, nb::rv_policy::take_ownership));
      },
      "mol"_a, "confId"_a = -1,
      "keepOverlapAndHamiltonianMatrices"_a = false,
      R"DOC(Runs an extended Hueckel calculation for a molecule.
The molecule should have at least one conformation

ARGUMENTS:
   - mol: molecule to use
   - confId: (optional) conformation to use
   - keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap
     and hamiltonian matrices in the EHTResults object

RETURNS: a 2-tuple:
   - a boolean indicating whether or not the calculation succeeded
   - an EHTResults object with the results
)DOC");
}