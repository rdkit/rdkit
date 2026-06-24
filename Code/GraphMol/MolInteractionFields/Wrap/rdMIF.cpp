//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/Wrap.h>
#include <RDBoost/PySequenceHolder.h>
#include <boost/python.hpp>
#include <ForceField/MMFF/Nonbonded.h>
#include <GraphMol/MolInteractionFields/MIFDescriptors.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/point.h>

#include <RDBoost/boost_numpy.h>

namespace python = boost::python;
using namespace RDMIF;

void wrap_mif();

BOOST_PYTHON_MODULE(rdMIF) {
  python::scope().attr("__doc__") =
      "Module containing functions for calculating molecular interaction fields (MIFs)\n\
  NOTE: This functionality is experimental and the API and/or results may change in future releases.";
  python::register_exception_translator<IndexErrorException>(
      &translate_index_error);
  python::register_exception_translator<ValueErrorException>(
      &translate_value_error);

  wrap_mif();
}

namespace RDMIF {

RDGeom::UniformRealValueGrid3D *constructGridHelper(const RDKit::ROMol &mol,
                                                    int confId, double margin,
                                                    double spacing) {
  return constructGrid(mol, confId, margin, spacing).release();
}

std::pair<std::vector<double>, std::vector<RDGeom::Point3D>>
extractChargesAndPositions(const python::object &charges,
                           const python::object &positions) {
  const auto pyPos = positions.ptr();
  const auto pyCharges = charges.ptr();
  if (!pyPos || !PySequence_Check(pyPos)) {
    throw_value_error("positions argument must be a sequence");
  }
  if (!pyCharges || !PySequence_Check(pyCharges)) {
    throw_value_error("charges argument must be a sequence");
  }
  auto nrows = PySequence_Size(pyPos);
  if (nrows != PySequence_Size(pyCharges)) {
    throw_value_error("positions and charges must have the same length");
  }

  auto extract_double = [](PyObject *obj, size_t i) {
    const auto dblObj = PySequence_GetItem(obj, i);
    double value = python::extract<double>(dblObj);
    Py_DecRef(dblObj);
    return value;
  };

  std::vector<RDGeom::Point3D> pos(nrows);
  std::vector<double> ch(nrows);
  for (unsigned int i = 0; i < nrows; ++i) {
    const auto pyXyz = PySequence_GetItem(pyPos, i);
    if (!pyXyz || !PySequence_Check(pyXyz) || PySequence_Size(pyXyz) != 3) {
      if (pyXyz) {
        Py_DecRef(pyXyz);
      }
      throw_value_error(
          "all elements in positions argument must be x,y,z sequences");
    }
    pos[i].x = extract_double(pyXyz, 0);
    pos[i].y = extract_double(pyXyz, 1);
    pos[i].z = extract_double(pyXyz, 2);
    ch[i] = extract_double(pyCharges, i);
    Py_DecRef(pyXyz);
  }
  return std::make_pair(std::move(ch), std::move(pos));
}

std::shared_ptr<Coulomb> makeAltCoulomb(const python::object &charges,
                                        const python::object &positions,
                                        double probecharge, bool absVal,
                                        double alpha, double cutoff) {
  const auto [ch, pos] = extractChargesAndPositions(charges, positions);
  return std::make_shared<Coulomb>(ch, pos, probecharge, absVal, alpha, cutoff);
}

std::shared_ptr<CoulombDielectric> makeAltCoulombDielectric(
    const python::object &charges, const python::object &positions,
    double probecharge, bool absVal, double alpha, double cutoff,
    double epsilon, double xi) {
  const auto [ch, pos] = extractChargesAndPositions(charges, positions);
  return std::make_shared<CoulombDielectric>(ch, pos, probecharge, absVal,
                                             alpha, cutoff, epsilon, xi);
}

python::tuple readCubeFile(const std::string &filename) {
  std::unique_ptr<RDGeom::UniformRealValueGrid3D> grd(
      new RDGeom::UniformRealValueGrid3D());
  auto res = readFromCubeFile(*grd, filename);
  boost::python::manage_new_object::apply<
      RDGeom::UniformRealValueGrid3D *>::type grdConverter;
  boost::python::manage_new_object::apply<RDKit::ROMol *>::type molConverter;
  return python::make_tuple(python::handle<>(grdConverter(grd.release())),
                            python::handle<>(molConverter(
                                static_cast<RDKit::ROMol *>(res.release()))));
}

struct mif_wrapper {
  static void wrap() {
    std::string docStringClass =
        "Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in\n\
        vacuum (no dielectric).\n\n";
    std::string docStringConst =
        "Constructor for Coulomb class.\n\n\
        ARGUMENTS:\n\
        - mol:           the molecule of interest\n\
        - confId:        the ID of the conformer to be used (defaults to -1)\n\
        - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
        - chargeKey      property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')\n\
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
        - cutoff         minimum cutoff distance [A] (defaults to 1.0)\n";
    std::string docStringConstAlt =
        "Alternative constructor for Coulomb class.\n\n\
        ARGUMENTS:\n\
        - charges:       array of partial charges of a molecule's atoms\n\
        - positions:     array of positions of a molecule's atoms\n\
        - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
        - cutoff         minimum cutoff distance [A] (defaults to 1.0)\n";
    std::string docString =
        "Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in\n\
        vacuum (no dielectric).\n\n\
        ARGUMENTS:\n\
        - x, y, z:   coordinates of probe position for energy calculation\n\
        - threshold: maximal distance until which interactions are calculated\n\
        RETURNS:\n\
        - electrostatic potential in [kJ mol^-1]\n";
    python::class_<Coulomb, std::shared_ptr<Coulomb>>(
        "Coulomb", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, double, bool,
                     const std::string &, double, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
             python::arg("chargeKey") = "_GasteigerCharge",
             python::arg("softcoreParam") = 0.0, python::arg("cutoff") = 1.0),
            docStringConst.c_str()))
        .def("__init__",
             python::make_constructor(
                 makeAltCoulomb, python::default_call_policies(),
                 (python::arg("charges"), python::arg("positions"),
                  python::arg("probeCharge") = 1.0,
                  python::arg("absVal") = false,
                  python::arg("softcoreParam") = 0.0,
                  python::arg("cutoff") = 1.0)),
             docStringConstAlt.c_str())
        .def("__call__", &Coulomb::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringClass =
        "Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in\n\
        by taking a distance-dependent dielectric into account.\n\
        Same energy term as used in GRID MIFs.\n\
        References:\n\
        - J. Med. Chem. 1985, 28, 849.\n\
        - J. Comp. Chem. 1983, 4, 187.\n\n";
    docStringConst =
        "Constructor for CoulombDielectric class.\n\n\
        ARGUMENTS:\n\
        - mol:           the molecule of interest\n\
        - confId:        the ID of the conformer to be used (defaults to -1)\n\
        - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
        - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
        - chargeKey       property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')\n\
        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
        - cutoff         minimum cutoff distance [A] (defaults to 1.0)\n\
        - epsilon        relative permittivity of solvent (defaults to 80.0)\n\
        - xi             relative permittivity of solute (defaults to 4.0)\n";
    docStringConstAlt =
        "Alternative constructor for CoulombDielectric class.\n\n\
        ARGUMENTS:\n\
      - charges:       array of partial charges of a molecule's atoms\n\
      - positions:     array of positions of a molecule's atoms\n\
      - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
      - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
      - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
      - cutoff         minimum cutoff distance [A] (defaults to 1.0)\n\
      - epsilon        relative permittivity of solvent (defaults to 80.0)\n\
      - xi             relative permittivity of solute (defaults to 4.0)\n";
    docString =
        "Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in\n\
        by taking a distance-dependent dielectric into account.\n\n\
        ARGUMENTS:\n\
        - x, y, z:   coordinates of probe position for energy calculation\n\
        - threshold: maximal distance until which interactions are calculated\n\
        RETURNS:\n\
        - electrostatic potential in [kJ mol^-1]\n";
    python::class_<CoulombDielectric, std::shared_ptr<CoulombDielectric>>(
        "CoulombDielectric", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, double, bool,
                     const std::string &, double, double, double, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
             python::arg("chargeKey") = "_GasteigerCharge",
             python::arg("softcoreParam") = 0.0, python::arg("cutoff") = 1.0,
             python::arg("epsilon") = 80.0, python::arg("xi") = 4.0),
            docStringConst.c_str()))
        .def("__init__",
             python::make_constructor(
                 makeAltCoulombDielectric, python::default_call_policies(),
                 (python::arg("charges"), python::arg("positions"),
                  python::arg("probeCharge") = 1.0,
                  python::arg("absVal") = false,
                  python::arg("softcoreParam") = 0.0,
                  python::arg("cutoff") = 1.0, python::arg("epsilon") = 80.0,
                  python::arg("xi") = 4.0)),
             docStringConstAlt.c_str())
        .def(python::init<const std::string &>(
            python::args("self", "pklString")))
        .def("__call__", &CoulombDielectric::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringClass =
        "Class for calculating van der Waals interactions between molecule and a probe at a gridpoint\
        based on the MMFF forcefield.\n";
    docStringConst =
        "ARGUMENTS:\n\
        - mol           molecule object\n\
        - confId        conformation id which is used to get positions of atoms (default=-1)\n\
        - probeAtomType MMFF94 atom type for the probe atom (default=6, sp3 oxygen)\n\
        - cutoff        minimum cutoff distance [A] (default:1.0)\n\
        - scaling       scaling of VdW parameters to take hydrogen bonds into account (default=False)\n";
    docString =
        "Calculates the van der Waals interaction between molecule and a probe at a gridpoint.\n\n\
        ARGUMENTS:\n\
        - x, y, z:   coordinates of probe position for energy calculation\n\
        - threshold: maximal distance until which interactions are calculated\n\
        RETURNS:\n\
        - van der Waals potential in [kJ mol^-1]\n";
    python::class_<MMFFVdWaals, std::shared_ptr<MMFFVdWaals>,
                   boost::noncopyable>(
        "MMFFVdWaals", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, unsigned int, bool, double>(
            (python::arg("self"), python::arg("mol"),
             python::arg("confId") = -1, python::arg("probeAtomType") = 6,
             python::arg("scaling") = false, python::arg("cutoff") = 1.0),
            docStringConst.c_str()))
        .def("__call__", &MMFFVdWaals::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringClass =
        "Class for calculating van der Waals interactions between molecule and a probe at a gridpoint\
        based on the UFF forcefield.\n";
    docStringConst =
        "ARGUMENTS:\n\
        - mol           molecule object\n\
        - confId        conformation id which is used to get positions of atoms (default=-1)\n\
        - probeAtomType UFF atom type for the probe atom (default='O_3', sp3 oxygen)\n\
        - cutoff        minimum cutoff distance [A] (default:1.0)\n";
    python::class_<UFFVdWaals, std::shared_ptr<UFFVdWaals>, boost::noncopyable>(
        "UFFVdWaals", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, const std::string &, double>(
            (python::arg("self"), python::arg("mol"),
             python::arg("confId") = -1, python::arg("probeAtomType") = "O_3",
             python::arg("cutoff") = 1.0),
            docStringConst.c_str()))
        .def("__call__", &UFFVdWaals::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringClass =
        "Class for calculation of hydrogen bonding energy between a probe and a molecule.\n\n\
        Similar to GRID hydrogen bonding descriptors.\n\
        References:\n\
        - J.Med.Chem. 1989, 32, 1083.\n\
        - J.Med.Chem. 1993, 36, 140.\n\
        - J.Med.Chem. 1993, 36, 148.\n";
    docStringConst =
        "Constructor for HBond class.\n\n\
        ARGUMENTS:\n\
        - mol:           the molecule of interest\n\
        - confId:        the ID of the conformer to be used (defaults to -1)\n\
        - probeAtomType: atom type for the probe atom (either 'OH', 'O', 'NH' or 'N') (defaults to 'OH')\n\
        - fixed:         for some groups, two different angle dependencies are defined:\n\
                         one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)\n\
                         into account and one for strictly fixed conformations\n\
                         if True, strictly fixed conformations (defaults to True)\n\
        - cutoff         minimum cutoff distance [A] (defaults to 1.0)\n";
    docString =
        "Calculates the hydrogen bonding energy between probe and molecule in\n\n\
        ARGUMENTS:\n\
        - x, y, z:   coordinates of probe position for energy calculation\n\
        - threshold: maximal distance until which interactions are calculated\n\
        RETURNS:\n\
        hydrogen bonding energy in [kJ mol^-1]\n";
    python::class_<HBond, std::shared_ptr<HBond>>(
        "HBond", docStringClass.c_str(),
        python::init<RDKit::ROMol &, int, const std::string &, bool, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeAtomType") = "OH", python::arg("fixed") = true,
             python::arg("cutoff") = 1.0),
            docStringConst.c_str()))
        .def("__call__", &HBond::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringClass =
        "Class for calculation of a hydrophilic potential of a molecule at a point.\n\n\
        The interaction energy of hydrogen and oxygen of water is calculated at each point as a \n\
        hydrogen bond interaction (either OH or O probe). The favored interaction is returned.\n";
    docStringConst =
        "Constructor for Hydrophilic class.\n\n\
        ARGUMENTS:\n\
        - mol:         the molecule of interest\n\
        - confId:      the ID of the conformer to be used (defaults to -1)\n\
        - fixed:       for some groups, two different angle dependencies are defined:\n\
                       one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)\n\
                       into account and one for strictly fixed conformations\n\
                       if True, strictly fixed conformations (defaults to True)\n\
        - cutoff       minimum cutoff distance [A] (default:1.0)\n";
    docString =
        "Calculates the hydrophilic field energy at a point.\n\n\
        ARGUMENTS:\n\
        - x, y, z:   coordinates of probe position for energy calculation\n\
        - threshold: maximal distance until which interactions are calculated\n\
        RETURNS:\n\
        hydrophilic field energy in [kJ mol^-1]\n";
    python::class_<Hydrophilic, std::shared_ptr<Hydrophilic>>(
        "Hydrophilic", docStringClass.c_str(),
        python::init<RDKit::ROMol &, int, bool, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("fixed") = true, python::arg("cutoff") = 1.0),
            docStringConst.c_str()))
        .def("__call__", &Hydrophilic::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docString =
        "Constructs a UniformRealValueGrid3D (3D grid with real values at gridpoints) fitting to a molecule.\n\n\
        ARGUMENTS:\n\
        - mol:     molecule of interest\n\
        - confId:  the ID of the conformer to be used (defaults to -1)\n\
        - margin:  minimum distance of molecule to surface of grid [A] (defaults to 5.0 A)\n\
        - spacing: grid spacing [A] (defaults to 0.5 A)\n";
    python::def("ConstructGrid", constructGridHelper,
                (python::arg("mol"), python::arg("confId") = -1,
                 python::arg("margin") = 5.0, python::arg("spacing") = 0.5),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

    docString =
        "Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.\n\n\
        ARGUMENTS:\n\
        - grid:      UniformRealValueGrid3D which get the MIF values\n\
        - descriptor:  Descriptor class which is used to calculate values\n";
    python::def("CalculateDescriptors", calculateDescriptors<Coulomb>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<CoulombDielectric>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<MMFFVdWaals>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<UFFVdWaals>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<HBond>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<Hydrophilic>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    docString =
        "Writes Grid to a file in Gaussian CUBE format.\n\n\
        ARGUMENTS:\n\
        - grid:      UniformRealValueGrid3D to be stored\n\
        - filename:  filename of file to be written\n\
        - mol:       associated molecule (defaults to None)\n\
        - confId:    the ID of the conformer to be used (defaults to -1)\n";
    python::def(
        "WriteToCubeFile", writeToCubeFile,
        (python::arg("grid"), python::arg("filename"),
         python::arg("mol") = python::object(), python::arg("confId") = -1),
        docString.c_str());

    docString =
        "Reads Grid from a file in Gaussian CUBE format.\n\n\
        ARGUMENTS:\n\
        - filename:  filename of file to be read\n\
        RETURNS:\n\
        a tuple where the first element is the grid and\n\
        the second element is the molecule object associated to the grid\n\
        (only atoms and coordinates, no bonds;\n\
        None if no molecule was associated to the grid)\n";
    python::def("ReadFromCubeFile", readCubeFile, (python::arg("filename")),
                docString.c_str());
  }
};
}  // namespace RDMIF

void wrap_mif() { mif_wrapper::wrap(); }
