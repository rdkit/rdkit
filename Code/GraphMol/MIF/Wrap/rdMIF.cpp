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
#include <GraphMol/MIF/MIFDescriptors.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/point.h>

#include <RDBoost/boost_numpy.h>

namespace python = boost::python;
using namespace RDMIF;

void wrap_mif();

BOOST_PYTHON_MODULE(rdMIF) {
  python::scope().attr("__doc__") =
      "Module containing functions for calculating molecular interaction fields (MIFs)";
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

Coulomb *make_coulomb(const python::object &charges,
                      const python::object &positions, double probecharge,
                      bool absVal, double alpha, double cutoff) {
  PyObject *pyObj = positions.ptr();
  std::vector<RDGeom::Point3D> pos;
  std::vector<double> ch;
  unsigned int nrows;
  if (PySequence_Check(pyObj)) {
    nrows = PySequence_Size(pyObj);
    if (nrows <= 0) throw_value_error("Empty sequence passed in");
    python::extract<RDGeom::Point3D> ptOk(positions[0]);
    if (!ptOk.check()) {
      for (unsigned int i = 0; i < nrows; i++) {
        PySequenceHolder<double> row(positions[i]);
        if (row.size() != 3)
          throw_value_error("Wrong number of entries in the list of lists");
        RDGeom::Point3D pt(row[0], row[1], row[2]);
        pos.push_back(pt);
        ch.push_back(python::extract<double>(charges[i]));
      }
    } else {
      for (unsigned int i = 0; i < nrows; i++) {
        python::extract<RDGeom::Point3D> pt(positions[i]);
        ch.push_back(python::extract<double>(charges[i]));
        if (pt.check()) {
          pos.push_back(pt);
        } else {
          throw_value_error("non-Point3D found in sequence of points");
        }
      }
    }
  }

  Coulomb *coul = new Coulomb(ch, pos, probecharge, absVal, alpha, cutoff);
  return coul;
}

CoulombDielectric *make_coulomb_dielectric(const python::object &charges,
                                           const python::object &positions,
                                           double probecharge, bool absVal,
                                           double alpha, double cutoff,
                                           double epsilon, double xi) {
  PyObject *pyObj = positions.ptr();
  std::vector<RDGeom::Point3D> pos;
  std::vector<double> ch;
  unsigned int nrows;
  if (PySequence_Check(pyObj)) {
    nrows = PySequence_Size(pyObj);
    if (nrows <= 0) throw_value_error("Empty sequence passed in");
    python::extract<RDGeom::Point3D> ptOk(positions[0]);
    if (!ptOk.check()) {
      for (unsigned int i = 0; i < nrows; i++) {
        PySequenceHolder<double> row(positions[i]);
        if (row.size() != 3)
          throw_value_error("Wrong number of entries in the list of lists");
        RDGeom::Point3D pt(row[0], row[1], row[2]);
        pos.push_back(pt);
        ch.push_back(python::extract<double>(charges[i]));
      }
    } else {
      for (unsigned int i = 0; i < nrows; i++) {
        python::extract<RDGeom::Point3D> pt(positions[i]);
        ch.push_back(python::extract<double>(charges[i]));
        if (pt.check()) {
          pos.push_back(pt);
        } else {
          throw_value_error("non-Point3D found in sequence of points");
        }
      }
    }
  }

  CoulombDielectric *coul = new CoulombDielectric(ch, pos, probecharge, absVal,
                                                  alpha, cutoff, epsilon, xi);
  return coul;
}

VdWaals *constructVdWMMFF(RDKit::ROMol &mol, int confId,
                          unsigned int probeAtomType, bool scaling,
                          double cutoff) {
  VdWaals *res =
      new VdWaals(mol, confId, probeAtomType, "", "MMFF94", scaling, cutoff);
  return res;
}

VdWaals *constructVdWUFF(RDKit::ROMol &mol, int confId,
                         const std::string &probeAtomType, double cutoff) {
  VdWaals *res =
      new VdWaals(mol, confId, 0, probeAtomType, "UFF", false, cutoff);
  return res;
}

RDKit::ROMol *readCubeFile(RDGeom::UniformRealValueGrid3D &grd,
                           const std::string &filename) {
  auto res = readFromCubeFile(grd, filename);
  return res.release();
}

struct mif_wrapper {
  static void wrap() {
    std::string docStringClass =
        "Class for calculation of the closest distance of a point to a molecule.\n\n";
    std::string docStringConst =
        "Constructor for DistaceToClosestAtom class.\n\n\
				ARGUMENTS:\n\
				    - mol:    the molecule of interest\n\
				    - confId: the ID of the conformer to be used (defaults to -1)\n";
    std::string docString =
        "Calculates the closest distance from a point to a the molecule\n\n\
				ARGUMENTS:\n\
				    - pt: Point3D from which the distance to the molecule is calculated\n\
				RETURNS:\n\
				    - closest distance in [A]\n";
    python::class_<DistanceToClosestAtom>(
        "DistanceToClosestAtom", docStringClass.c_str(),
        python::init<const RDKit::ROMol, int>(
            (python::arg("mol"), python::arg("confId") = -1),
            docStringConst.c_str()))
        .def("__call__", &DistanceToClosestAtom::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z")),
             docString.c_str());

    docStringClass =
        "Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in\n\
						vaccuum (no dielectric).\n\n";
    docStringConst =
        "Constructor for Coulomb class.\n\n\
				ARGUMENTS:\n\
				- mol:           the molecule of interest\n\
				- confId:        the ID of the conformer to be used (defaults to -1)\n\
				- probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
			      	- absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
			        - chargeKey	 property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')\n\
			        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
			        - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n";
    docString =
        "Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in\n\
						vaccuum (no dielectric).\n\n\
				ARGUMENTS:\n\
				    - x, y, z:	 coordinates of probe position for energy calculation\n\
                                    - threshold: maximal distance until which interactions are calculated\n\
				RETURNS:\n\
					- electrostatic potential in [kJ mol^-1]\n";
    python::class_<Coulomb>(
        "Coulomb", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, double, bool,
                     const std::string &, double, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
             python::arg("chargeKey") = "_GasteigerCharge",
             python::arg("softcoreParam") = 0.0,
             python::arg("cutoffDist") = 1.0),
            docStringConst.c_str()))
        .def("__call__", &Coulomb::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringConst =
        "Alternative constructor for Coulomb class.\n\n\
						ARGUMENTS:\n\
						- charges:       array of partial charges of a molecule's atoms\n\
						- positions:     array of positions of a molecule's atoms\n\
					        - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
					      	- absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
					        - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
					        - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n";
    python::def(
        "Coulomb_", make_coulomb,
        (python::arg("charges"), python::arg("positions"),
         python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
         python::arg("softcoreParam") = 0.0, python::arg("cutoffDist") = 1.0),
        docStringConst.c_str(),
        python::return_value_policy<python::manage_new_object>());

    docStringClass =
        "Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in\n\
						by taking a distance-dependent dielectric into account:\n\
                        Same energy term as used in GRID MIFs\n\
			            References:\n\
                        J. Med. Chem. 1985, 28, 849.\n\
                        J. Comp. Chem. 1983, 4, 187.\n\n";
    docStringConst =
        "Constructor for CoulombDielectric class.\n\n\
				ARGUMENTS:\n\
				    - mol:           the molecule of interest\n\
				    - confId:        the ID of the conformer to be used (defaults to -1)\n\
				    - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
			      	    - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
			            - chargeKey	     property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')\n\
			            - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
			            - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n\
				    - epsilon        relative permittivity of solvent (defaults to 80.0)\n\
	  	                    - xi             relative permittivity of solute (defaults to 4.0)\n";
    docString =
        "Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in\n\
						by taking a distance-dependent dielectric into account.\n\n\
				ARGUMENTS:\n\
				    - x, y, z:	 coordinates of probe position for energy calculation\n\
                                    - threshold: maximal distance until which interactions are calculated\n\
                                RETURNS:\n\
					- electrostatic potential in [kJ mol^-1]\n";
    python::class_<CoulombDielectric>(
        "CoulombDielectric", docStringClass.c_str(),
        python::init<const RDKit::ROMol &, int, double, bool,
                     const std::string &, double, double, double, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
             python::arg("chargeKey") = "_GasteigerCharge",
             python::arg("softcoreParam") = 0.0,
             python::arg("cutoffDist") = 1.0, python::arg("epsilon") = 80.0,
             python::arg("xi") = 4.0),
            docStringConst.c_str()))
        .def("__call__", &CoulombDielectric::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringConst =
        "Alternative constructor for CoulombDielectric class.\n\n\
								ARGUMENTS:\n\
								    - charges:       array of partial charges of a molecule's atoms\n\
								    - positions:     array of positions of a molecule's atoms\n\
							            - probeCharge    charge of probe [e] (defaults to 1.0 e)\n\
							      	    - absVal:        if True, absolute values of interactions are calculated (defaults to False)\n\
							            - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)\n\
							            - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n\
								    - epsilon 	     relative permittivity of solvent (defaults to 80.0)\n\
								    - xi             relative permittivity of solute (defaults to 4.0)\n";
    python::def(
        "CoulombDielectric_", make_coulomb_dielectric,
        (python::arg("charges"), python::arg("positions"),
         python::arg("probeCharge") = 1.0, python::arg("absVal") = false,
         python::arg("softcoreParam") = 0.0, python::arg("cutoffDist") = 1.0,
         python::arg("epsilon") = 80.0, python::arg("xi") = 4.0),
        docStringConst.c_str(),
        python::return_value_policy<python::manage_new_object>());

    docStringClass =
        "Class for calculation van der Waals interaction between molecule and a probe at a gridpoint.\n\n";
    docString =
        "Calculates the van der Waals interaction between molecule and a probe at a gridpoint.\n\n\
				ARGUMENTS:\n\
				    - x, y, z:	 coordinates of probe position for energy calculation\n\
                                    - threshold: maximal distance until which interactions are calculated\n\
                                RETURNS:\n\
					- van der Waals potential in [kJ mol^-1]\n";
    python::class_<VdWaals>("VdWaals", "",
                            python::init<>("Default Constructor"))
        .def("__call__", &VdWaals::operator(),
             (python::arg("x"), python::arg("y"), python::arg("z"),
              python::arg("threshold")),
             docString.c_str());

    docStringConst =
        "Constructs VdWaals class which uses MMFF94 force field parameters.\n\n\
				ARGUMENTS:\n\
				    - mol:           the molecule of interest\n\
				    - confId:        the ID of the conformer to be used (defaults to -1)\n\
				    - probeType      MMFF94 atom type (integer) used as probe atom (defaults to 6)\n\
			      	    - scaling:       scales interaction to take hydrogen bonds into account (MMFF94-specific) (defaults to False)\n\
			            - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n";
    python::def("ConstructVdWaalsMMFF", &constructVdWMMFF,
                (python::arg("mol"), python::arg("confId") = -1,
                 python::arg("probeType") = 6, python::arg("scaling") = false,
                 python::arg("cutoffDist") = 1.0),
                docStringConst.c_str(),
                python::return_value_policy<python::manage_new_object>());

    docStringConst =
        "Constructs VdWaals class which uses UFF force field parameters.\n\n\
				ARGUMENTS:\n\
				    - mol:           the molecule of interest\n\
				    - confId:        the ID of the conformer to be used (defaults to -1)\n\
				    - probeType	     UFF atom type (string) used as probe atom (defaults to 'O_3')\n\
			            - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n";
    python::def(
        "ConstructVdWaalsUFF", &constructVdWUFF,
        (python::arg("mol"), python::arg("confId") = -1,
         python::arg("probeType") = "O_3", python::arg("cutoffDist") = 1.0),
        docStringConst.c_str(),
        python::return_value_policy<python::manage_new_object>());

    docStringClass =
        "Class for calculation of hydrogen bonding energy between a probe and a molecule.\n\n\
				Similar to GRID hydrogen bonding descriptors\n\
				References:\n\
				  - J.Med.Chem. 1989, 32, 1083.\n\
				  - J.Med.Chem. 1993, 36, 140.\n\
				  - J.Med.Chem. 1993, 36, 148.\n";
    docStringConst =
        "Constructor for HBond class.\n\n\
				ARGUMENTS:\n\
				    - mol:           the molecule of interest\n\
				    - confId:        the ID of the conformer to be used (defaults to -1)\n\
	                            - probeType:     atom type for the probe atom (either 'OH', 'O', 'NH' or 'N') (defaults to 'OH')\n\
	                            - fixed:         for some groups, two different angle dependencies are defined:\n\
	      	                         one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)\n\
	      	                         into account and one for strictly fixed conformations\n\
	      	                         if True, strictly fixed conformations (defaults to True)\n\
			            - cutoffDist     minimum cutoff distance [A] (defaults to 1.0)\n";
    docString =
        "Calculates the hydrogen bonding energy between probe and molecule in\n\n\
				ARGUMENTS:\n\
				    - x, y, z:	 coordinates of probe position for energy calculation\n\
                                    - threshold: maximal distance until which interactions are calculated\n\
                                RETURNS:\n\
					hydrogen bonding energy in [kJ mol^-1]\n";
    python::class_<HBond>(
        "HBond", docStringClass.c_str(),
        python::init<RDKit::ROMol &, int, const std::string &, bool, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("probeType") = "OH", python::arg("fixed") = true,
             python::arg("cutoffDist") = 1.0),
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
				    - mol:               the molecule of interest\n\
				    - confId:        	 the ID of the conformer to be used (defaults to -1)\n\
	                            - fixed:         	 for some groups, two different angle dependencies are defined:\n\
	      	                             one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)\n\
	      	                             into account and one for strictly fixed conformations\n\
	      	                             if True, strictly fixed conformations (defaults to True)\n\
	                            - cutoffDist	 minimum cutoff distance [A] (default:1.0)\n";
    docString =
        "Calculates the hydrophilic field energy at a point.\n\n\
				ARGUMENTS:\n\
				    - x, y, z:	 coordinates of probe position for energy calculation\n\
                                    - threshold: maximal distance until which interactions are calculated\n\
                                RETURNS:\n\
					hydrophilic field energy in [kJ mol^-1]\n";
    python::class_<Hydrophilic>(
        "Hydrophilic", docStringClass.c_str(),
        python::init<RDKit::ROMol &, int, bool, double>(
            (python::arg("mol"), python::arg("confId") = -1,
             python::arg("fixed") = true, python::arg("cutoffDist") = 1.0),
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
				    - grid: 	   UniformRealValueGrid3D which get the MIF values\n\
				    - descriptor:  Descriptor class which is used to calculate values\n";
    python::def(
        "CalculateDescriptors", calculateDescriptors<DistanceToClosestAtom>,
        (python::arg("grid"), python::arg("descriptor")), docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<Coulomb>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<CoulombDielectric>,
                (python::arg("grid"), python::arg("descriptor"),
                 python::arg("threshold") = -1.0),
                docString.c_str());

    python::def("CalculateDescriptors", calculateDescriptors<VdWaals>,
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
				    - grid: 	UniformRealValueGrid3D to be stored\n\
				    - mol:	respective molecule\n\
				    - filename:	filename of file to be written\n\
				    - confId:	the ID of the conformer to be used (defaults to -1)\n";
    python::def("WriteToCubeFile", writeToCubeFile,
                (python::arg("grid"), python::arg("mol"),
                 python::arg("filename"), python::arg("confId") = -1),
                docString.c_str());

    docString =
        "Reads Grid from a file in Gaussian CUBE format.\n\n\
				ARGUMENTS:\n\
				    - grid: 	UniformRealValueGrid3D where data is read in\n\
				 	- filename:	filename of file to be read\n\
				RETURNS:\n\
				    a molecule object (only atoms and coordinates, no bonds!)\n";
    python::def("ReadFromCubeFile", readCubeFile,
                (python::arg("grid"), python::arg("filename")),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());
  }
};
}  // namespace RDMIF

void wrap_mif() { mif_wrapper::wrap(); }
