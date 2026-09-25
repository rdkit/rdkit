//
//  Copyright (c) 2014-2026, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <GraphMol/MolInteractionFields/MIFDescriptors.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/point.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDMIF;

namespace {

std::pair<std::vector<double>, std::vector<RDGeom::Point3D>>
extractChargesAndPositions(const nb::object &charges,
                           const nb::object &positions) {
  if (positions.is_none() || !nb::isinstance<nb::sequence>(positions)) {
    throw ValueErrorException("positions argument must be a sequence");
  }
  if (charges.is_none() || !nb::isinstance<nb::sequence>(charges)) {
    throw ValueErrorException("charges argument must be a sequence");
  }
  nb::sequence pyPos = nb::cast<nb::sequence>(positions);
  nb::sequence pyCharges = nb::cast<nb::sequence>(charges);
  auto nrows = nb::len(pyPos);
  if (nrows != nb::len(pyCharges)) {
    throw ValueErrorException("positions and charges must have the same length");
  }

  std::vector<RDGeom::Point3D> pos(nrows);
  std::vector<double> ch(nrows);
  for (size_t i = 0; i < nrows; ++i) {
    nb::object pyXyz = pyPos[i];
    if (!nb::isinstance<nb::sequence>(pyXyz) || nb::len(pyXyz) != 3) {
      throw ValueErrorException(
          "all elements in positions argument must be x,y,z sequences");
    }
    nb::sequence xyz = nb::cast<nb::sequence>(pyXyz);
    pos[i].x = nb::cast<double>(xyz[0]);
    pos[i].y = nb::cast<double>(xyz[1]);
    pos[i].z = nb::cast<double>(xyz[2]);
    ch[i] = nb::cast<double>(pyCharges[i]);
  }
  return std::make_pair(std::move(ch), std::move(pos));
}

std::tuple<RDGeom::UniformRealValueGrid3D *, RDKit::ROMol *> readCubeFileHelper(
    const std::string &filename) {
  std::unique_ptr<RDGeom::UniformRealValueGrid3D> grd(
      new RDGeom::UniformRealValueGrid3D());
  auto res = readFromCubeFile(*grd, filename);
  return std::make_tuple(grd.release(),
                         static_cast<RDKit::ROMol *>(res.release()));
}

}  // namespace

NB_MODULE(rdMIF, m) {
  m.doc() =
      R"DOC(Module containing functions for calculating molecular interaction fields (MIFs)
NOTE: This functionality is experimental and the API and/or results may change in future releases.)DOC";

  nb::exception<ValueErrorException>(m, "MIFValueError", PyExc_ValueError);
  nb::exception<IndexErrorException>(m, "MIFIndexError", PyExc_IndexError);

  nb::class_<Coulomb>(m, "Coulomb",
                      R"DOC(Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
        vacuum (no dielectric).
)DOC")
      .def(nb::init<const RDKit::ROMol &, int, double, bool,
                    const std::string &, double, double>(),
           "mol"_a, "confId"_a = -1, "probeCharge"_a = 1.0,
           "absVal"_a = false, "chargeKey"_a = "_GasteigerCharge",
           "softcoreParam"_a = 0.0, "cutoff"_a = 1.0,
           R"DOC(Constructor for Coulomb class.

        ARGUMENTS:
        - mol:           the molecule of interest
        - confId:        the ID of the conformer to be used (defaults to -1)
        - probeCharge    charge of probe [e] (defaults to 1.0 e)
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)
        - chargeKey      property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
        - cutoff         minimum cutoff distance [A] (defaults to 1.0))DOC")
      .def("__init__",
           [](Coulomb &self, const nb::object &charges,
              const nb::object &positions, double probeCharge, bool absVal,
              double softcoreParam, double cutoff) {
             const auto [ch, pos] =
                 extractChargesAndPositions(charges, positions);
             new (&self) Coulomb(ch, pos, probeCharge, absVal, softcoreParam,
                                 cutoff);
           },
           "charges"_a, "positions"_a, "probeCharge"_a = 1.0,
           "absVal"_a = false, "softcoreParam"_a = 0.0, "cutoff"_a = 1.0,
           R"DOC(Alternative constructor for Coulomb class.

        ARGUMENTS:
        - charges:       array of partial charges of a molecule's atoms
        - positions:     array of positions of a molecule's atoms
        - probeCharge    charge of probe [e] (defaults to 1.0 e)
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
        - cutoff         minimum cutoff distance [A] (defaults to 1.0))DOC")
      .def("__call__", &Coulomb::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
        vacuum (no dielectric).

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        - electrostatic potential in [kJ mol^-1])DOC");

  nb::class_<CoulombDielectric>(
      m, "CoulombDielectric",
      R"DOC(Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
        by taking a distance-dependent dielectric into account.
        Same energy term as used in GRID MIFs.
        References:
        - J. Med. Chem. 1985, 28, 849.
        - J. Comp. Chem. 1983, 4, 187.
)DOC")
      .def(nb::init<const RDKit::ROMol &, int, double, bool,
                    const std::string &, double, double, double, double>(),
           "mol"_a, "confId"_a = -1, "probeCharge"_a = 1.0,
           "absVal"_a = false, "chargeKey"_a = "_GasteigerCharge",
           "softcoreParam"_a = 0.0, "cutoff"_a = 1.0, "epsilon"_a = 80.0,
           "xi"_a = 4.0,
           R"DOC(Constructor for CoulombDielectric class.

        ARGUMENTS:
        - mol:           the molecule of interest
        - confId:        the ID of the conformer to be used (defaults to -1)
        - probeCharge    charge of probe [e] (defaults to 1.0 e)
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)
        - chargeKey       property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
        - cutoff         minimum cutoff distance [A] (defaults to 1.0)
        - epsilon        relative permittivity of solvent (defaults to 80.0)
        - xi             relative permittivity of solute (defaults to 4.0))DOC")
      .def("__init__",
           [](CoulombDielectric &self, const nb::object &charges,
              const nb::object &positions, double probeCharge, bool absVal,
              double softcoreParam, double cutoff, double epsilon, double xi) {
             const auto [ch, pos] =
                 extractChargesAndPositions(charges, positions);
             new (&self) CoulombDielectric(ch, pos, probeCharge, absVal,
                                           softcoreParam, cutoff, epsilon, xi);
           },
           "charges"_a, "positions"_a, "probeCharge"_a = 1.0,
           "absVal"_a = false, "softcoreParam"_a = 0.0, "cutoff"_a = 1.0,
           "epsilon"_a = 80.0, "xi"_a = 4.0,
           R"DOC(Alternative constructor for CoulombDielectric class.

      - charges:       array of partial charges of a molecule's atoms
      - positions:     array of positions of a molecule's atoms
      - probeCharge    charge of probe [e] (defaults to 1.0 e)
      - absVal:        if True, absolute values of interactions are calculated (defaults to False)
      - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
      - cutoff         minimum cutoff distance [A] (defaults to 1.0)
      - epsilon        relative permittivity of solvent (defaults to 80.0)
      - xi             relative permittivity of solute (defaults to 4.0))DOC")
      .def("__call__", &CoulombDielectric::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
        by taking a distance-dependent dielectric into account.

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        - electrostatic potential in [kJ mol^-1])DOC");

  nb::class_<MMFFVdWaals>(
      m, "MMFFVdWaals",
      R"DOC(Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the MMFF forcefield.
)DOC")
      .def(nb::init<const RDKit::ROMol &, int, unsigned int, bool, double>(),
           "mol"_a, "confId"_a = -1, "probeAtomType"_a = 6,
           "scaling"_a = false, "cutoff"_a = 1.0,
           R"DOC(ARGUMENTS:
        - mol           molecule object
        - confId        conformation id which is used to get positions of atoms (default=-1)
        - probeAtomType MMFF94 atom type for the probe atom (default=6, sp3 oxygen)
        - cutoff        minimum cutoff distance [A] (default:1.0)
        - scaling       scaling of VdW parameters to take hydrogen bonds into account (default=False))DOC")
      .def("__call__", &MMFFVdWaals::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the van der Waals interaction between molecule and a probe at a gridpoint.

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        - van der Waals potential in [kJ mol^-1])DOC");

  nb::class_<UFFVdWaals>(
      m, "UFFVdWaals",
      R"DOC(Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the UFF forcefield.
)DOC")
      .def(nb::init<const RDKit::ROMol &, int, const std::string &, double>(),
           "mol"_a, "confId"_a = -1, "probeAtomType"_a = "O_3",
           "cutoff"_a = 1.0,
           R"DOC(ARGUMENTS:
        - mol           molecule object
        - confId        conformation id which is used to get positions of atoms (default=-1)
        - probeAtomType UFF atom type for the probe atom (default='O_3', sp3 oxygen)
        - cutoff        minimum cutoff distance [A] (default:1.0))DOC")
      .def("__call__", &UFFVdWaals::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the van der Waals interaction between molecule and a probe at a gridpoint.

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        - van der Waals potential in [kJ mol^-1])DOC");

  nb::class_<HBond>(m, "HBond",
                    R"DOC(Class for calculation of hydrogen bonding energy between a probe and a molecule.

        Similar to GRID hydrogen bonding descriptors.
        References:
        - J.Med.Chem. 1989, 32, 1083.
        - J.Med.Chem. 1993, 36, 140.
        - J.Med.Chem. 1993, 36, 148.
)DOC")
      .def(nb::init<RDKit::ROMol &, int, const std::string &, bool, double>(),
           "mol"_a, "confId"_a = -1, "probeAtomType"_a = "OH",
           "fixed"_a = true, "cutoff"_a = 1.0,
           R"DOC(Constructor for HBond class.

        ARGUMENTS:
        - mol:           the molecule of interest
        - confId:        the ID of the conformer to be used (defaults to -1)
        - probeAtomType: atom type for the probe atom (either 'OH', 'O', 'NH' or 'N') (defaults to 'OH')
        - fixed:         for some groups, two different angle dependencies are defined:
                         one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)
                         into account and one for strictly fixed conformations
                         if True, strictly fixed conformations (defaults to True)
        - cutoff         minimum cutoff distance [A] (defaults to 1.0))DOC")
      .def("__call__", &HBond::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the hydrogen bonding energy between probe and molecule in

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        hydrogen bonding energy in [kJ mol^-1])DOC");

  nb::class_<Hydrophilic>(
      m, "Hydrophilic",
      R"DOC(Class for calculation of a hydrophilic potential of a molecule at a point.

        The interaction energy of hydrogen and oxygen of water is calculated at each point as a
        hydrogen bond interaction (either OH or O probe). The favored interaction is returned.
)DOC")
      .def(nb::init<RDKit::ROMol &, int, bool, double>(), "mol"_a,
           "confId"_a = -1, "fixed"_a = true, "cutoff"_a = 1.0,
           R"DOC(Constructor for Hydrophilic class.

        ARGUMENTS:
        - mol:         the molecule of interest
        - confId:      the ID of the conformer to be used (defaults to -1)
        - fixed:       for some groups, two different angle dependencies are defined:
                       one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)
                       into account and one for strictly fixed conformations
                       if True, strictly fixed conformations (defaults to True)
        - cutoff       minimum cutoff distance [A] (default:1.0))DOC")
      .def("__call__", &Hydrophilic::operator(), "x"_a, "y"_a, "z"_a,
           "threshold"_a,
           R"DOC(Calculates the hydrophilic field energy at a point.

        ARGUMENTS:
        - x, y, z:   coordinates of probe position for energy calculation
        - threshold: maximal distance until which interactions are calculated
        RETURNS:
        hydrophilic field energy in [kJ mol^-1])DOC");

  m.def(
      "ConstructGrid",
      [](const RDKit::ROMol &mol, int confId, double margin, double spacing) {
        return constructGrid(mol, confId, margin, spacing).release();
      },
      "mol"_a, "confId"_a = -1, "margin"_a = 5.0, "spacing"_a = 0.5,
      R"DOC(Constructs a UniformRealValueGrid3D (3D grid with real values at gridpoints) fitting to a molecule.

        ARGUMENTS:
        - mol:     molecule of interest
        - confId:  the ID of the conformer to be used (defaults to -1)
        - margin:  minimum distance of molecule to surface of grid [A] (defaults to 5.0 A)
        - spacing: grid spacing [A] (defaults to 0.5 A))DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid, const Coulomb &descriptor,
         double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid,
         const CoulombDielectric &descriptor, double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid, const MMFFVdWaals &descriptor,
         double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid, const UFFVdWaals &descriptor,
         double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid, const HBond &descriptor,
         double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "CalculateDescriptors",
      [](RDGeom::UniformRealValueGrid3D &grid, const Hydrophilic &descriptor,
         double threshold) {
        calculateDescriptors(grid, descriptor, threshold);
      },
      "grid"_a, "descriptor"_a, "threshold"_a = -1.0,
      R"DOC(Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D which get the MIF values
        - descriptor:  Descriptor class which is used to calculate values)DOC");

  m.def(
      "WriteToCubeFile",
      [](const RDGeom::UniformRealValueGrid3D &grid, const std::string &filename,
         const nb::object &mol, int confId) {
        const RDKit::ROMol *molPtr = nullptr;
        if (!mol.is_none()) {
          molPtr = nb::cast<const RDKit::ROMol *>(mol);
        }
        writeToCubeFile(grid, filename, molPtr, confId);
      },
      "grid"_a, "filename"_a, "mol"_a = nb::none(), "confId"_a = -1,
      R"DOC(Writes Grid to a file in Gaussian CUBE format.

        ARGUMENTS:
        - grid:      UniformRealValueGrid3D to be stored
        - filename:  filename of file to be written
        - mol:       associated molecule (defaults to None)
        - confId:    the ID of the conformer to be used (defaults to -1))DOC");

  m.def(
      "ReadFromCubeFile", readCubeFileHelper, "filename"_a,
      R"DOC(Reads Grid from a file in Gaussian CUBE format.

        ARGUMENTS:
        - filename:  filename of file to be read
        RETURNS:
        a tuple where the first element is the grid and
        the second element is the molecule object associated to the grid
        (only atoms and coordinates, no bonds;
        None if no molecule was associated to the grid))DOC",
      nb::rv_policy::take_ownership);
}
