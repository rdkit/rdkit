//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>

#include <GraphMol/GraphMol.h>

#include <ForceField/ForceField.h>
#include <ForceField/nbWrap/PyForceField.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
int UFFHelper(ROMol &mol, int maxIters, double vdwThresh, int confId,
              bool ignoreInterfragInteractions) {
  nb::gil_scoped_release release;
  return UFF::UFFOptimizeMolecule(mol, maxIters, vdwThresh, confId,
                                  ignoreInterfragInteractions)
      .first;
}

nb::list UFFConfsHelper(ROMol &mol, int numThreads, int maxIters,
                        double vdwThresh, bool ignoreInterfragInteractions) {
  std::vector<std::pair<int, double>> res;
  {
    nb::gil_scoped_release release;
    UFF::UFFOptimizeMoleculeConfs(mol, res, numThreads, maxIters, vdwThresh,
                                  ignoreInterfragInteractions);
  }
  nb::list pyres;
  for (auto &itm : res) {
    pyres.append(nb::make_tuple(itm.first, itm.second));
  }
  return pyres;
}

nb::list MMFFConfsHelper(ROMol &mol, int numThreads, int maxIters,
                         std::string mmffVariant, double nonBondedThresh,
                         bool ignoreInterfragInteractions) {
  std::vector<std::pair<int, double>> res;
  {
    nb::gil_scoped_release release;
    MMFF::MMFFOptimizeMoleculeConfs(mol, res, numThreads, maxIters, mmffVariant,
                                    nonBondedThresh,
                                    ignoreInterfragInteractions);
  }
  nb::list pyres;
  for (auto &itm : res) {
    pyres.append(nb::make_tuple(itm.first, itm.second));
  }
  return pyres;
}

int FFHelper(ForceFields::PyForceField &ff, int maxIters) {
  nb::gil_scoped_release release;
  return ForceFieldsHelper::OptimizeMolecule(*ff.field, maxIters).first;
}

nb::list FFConfsHelper(ROMol &mol, ForceFields::PyForceField &ff,
                       int numThreads, int maxIters) {
  std::vector<std::pair<int, double>> res;
  {
    nb::gil_scoped_release release;
    ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff.field, res, numThreads,
                                             maxIters);
  }
  nb::list pyres;
  for (auto &itm : res) {
    pyres.append(nb::make_tuple(itm.first, itm.second));
  }
  return pyres;
}

std::unique_ptr<ForceFields::PyForceField> CreateEmptyForceFieldForMol(
    ROMol &mol, int confId = -1) {
  auto ff = ForceFieldsHelper::createEmptyForceFieldForMol(mol, confId);
  auto res = std::make_unique<ForceFields::PyForceField>(std::move(ff));
  res->initialize();
  return res;
}

std::unique_ptr<ForceFields::PyForceField> UFFGetMoleculeForceField(
    ROMol &mol, double vdwThresh = 10.0, int confId = -1,
    bool ignoreInterfragInteractions = true) {
  auto ff = std::unique_ptr<ForceFields::ForceField>(UFF::constructForceField(
      mol, vdwThresh, confId, ignoreInterfragInteractions));
  auto res = std::make_unique<ForceFields::PyForceField>(std::move(ff));
  res->initialize();
  return res;
}

bool UFFHasAllMoleculeParams(const ROMol &mol) {
  auto [types, foundAll] = UFF::getAtomTypes(mol);
  return foundAll;
}

int MMFFOptimizeMolecule(ROMol &mol, std::string mmffVariant = "MMFF94",
                         int maxIters = 200, double nonBondedThresh = 100.0,
                         int confId = -1,
                         bool ignoreInterfragInteractions = true) {
  int res = -1;

  MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
  if (mmffMolProperties.isValid()) {
    nb::gil_scoped_release release;
    std::unique_ptr<ForceFields::ForceField> ff(
        MMFF::constructForceField(mol, &mmffMolProperties, nonBondedThresh,
                                  confId, ignoreInterfragInteractions));
    ff->initialize();
    res = ff->minimize(maxIters);
  }

  return res;
}

unsigned int SanitizeMMFFMol(ROMol &mol) {
  return MMFF::sanitizeMMFFMol(static_cast<RWMol &>(mol));
};

std::unique_ptr<MMFF::MMFFMolProperties> GetMMFFMolProperties(
    ROMol &mol, std::string mmffVariant = "MMFF94",
    unsigned int mmffVerbosity = MMFF::MMFF_VERBOSITY_NONE) {
  auto mmffMolProperties = std::make_unique<MMFF::MMFFMolProperties>(
      mol, mmffVariant, mmffVerbosity);
  if (!mmffMolProperties->isValid()) {
    return nullptr;
  }
  return mmffMolProperties;
}

std::unique_ptr<ForceFields::PyForceField> MMFFGetMoleculeForceField(
    ROMol &mol, MMFF::MMFFMolProperties *MMFFMolProperties,
    double nonBondedThresh = 100.0, int confId = -1,
    bool ignoreInterfragInteractions = true) {
  if (!MMFFMolProperties) {
    return nullptr;
  }
  auto ff = std::unique_ptr<ForceFields::ForceField>(
      MMFF::constructForceField(mol, MMFFMolProperties, nonBondedThresh, confId,
                                ignoreInterfragInteractions));
  auto pyFF = std::make_unique<ForceFields::PyForceField>(std::move(ff));
  pyFF->initialize();
  return pyFF;
}

bool MMFFHasAllMoleculeParams(const ROMol &mol) {
  ROMol molCopy(mol);
  MMFF::MMFFMolProperties mmffMolProperties(molCopy);

  return mmffMolProperties.isValid();
}
};  // namespace RDKit

namespace ForceFields {
nb::object getUFFBondStretchParams(const RDKit::ROMol &mol,
                                   const unsigned int idx1,
                                   const unsigned int idx2) {
  ForceFields::UFF::UFFBond uffBondStretchParams;
  if (RDKit::UFF::getUFFBondStretchParams(mol, idx1, idx2,
                                          uffBondStretchParams)) {
    return nb::make_tuple(uffBondStretchParams.kb, uffBondStretchParams.r0);
  }
  return nb::none();
};

nb::object getUFFAngleBendParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3) {
  ForceFields::UFF::UFFAngle uffAngleBendParams;
  if (RDKit::UFF::getUFFAngleBendParams(mol, idx1, idx2, idx3,
                                        uffAngleBendParams)) {
    return nb::make_tuple(uffAngleBendParams.ka, uffAngleBendParams.theta0);
  }
  return nb::none();
};

nb::object getUFFTorsionParams(const RDKit::ROMol &mol, const unsigned int idx1,
                               const unsigned int idx2, const unsigned int idx3,
                               const unsigned int idx4) {
  ForceFields::UFF::UFFTor uffTorsionParams;
  if (RDKit::UFF::getUFFTorsionParams(mol, idx1, idx2, idx3, idx4,
                                      uffTorsionParams)) {
    return nb::float_(uffTorsionParams.V);
  }
  return nb::none();
};

nb::object getUFFInversionParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4) {
  ForceFields::UFF::UFFInv uffInversionParams;
  if (RDKit::UFF::getUFFInversionParams(mol, idx1, idx2, idx3, idx4,
                                        uffInversionParams)) {
    return nb::float_(uffInversionParams.K);
  }
  return nb::none();
};

nb::object getUFFVdWParams(const RDKit::ROMol &mol, const unsigned int idx1,
                           const unsigned int idx2) {
  ForceFields::UFF::UFFVdW uffVdWParams;
  if (RDKit::UFF::getUFFVdWParams(mol, idx1, idx2, uffVdWParams)) {
    return nb::make_tuple(uffVdWParams.x_ij, uffVdWParams.D_ij);
  }
  return nb::none();
};
}  // namespace ForceFields

NB_MODULE(rdForceFieldHelpers, m) {
  nb::module_::import_("rdkit.ForceField.rdForceField");

  m.doc() = "Module containing functions to handle force fields";

  m.def("UFFOptimizeMolecule", RDKit::UFFHelper, "mol"_a, "maxIters"_a = 200,
        "vdwThresh"_a = 10.0, "confId"_a = -1,
        "ignoreInterfragInteractions"_a = true,
        R"DOC(uses UFF to optimize a molecule's structure

ARGUMENTS:

 - mol : the molecule of interest
 - maxIters : the maximum number of iterations (defaults to 200)
 - vdwThresh : used to exclude long-range van der Waals interactions
               (defaults to 10.0)
 - confId : indicates which conformer to optimize
 - ignoreInterfragInteractions : if true, nonbonded terms between
               fragments will not be added to the forcefield.

RETURNS: 0 if the optimization converged, 1 if more iterations are required.
)DOC");

  m.def("UFFOptimizeMoleculeConfs", RDKit::UFFConfsHelper, "mol"_a,
        "numThreads"_a = 1, "maxIters"_a = 200, "vdwThresh"_a = 10.0,
        "ignoreInterfragInteractions"_a = true,
        R"DOC(uses UFF to optimize all of a molecule's conformations

ARGUMENTS:

 - mol : the molecule of interest
 - numThreads : the number of threads to use, only has an effect if the RDKit
                was built with thread support (defaults to 1)
                If set to zero, the max supported by the system will be used.
 - maxIters : the maximum number of iterations (defaults to 200)
 - vdwThresh : used to exclude long-range van der Waals interactions
               (defaults to 10.0)
 - ignoreInterfragInteractions : if true, nonbonded terms between
               fragments will not be added to the forcefield.

RETURNS: a list of (not_converged, energy) 2-tuples.
    If not_converged is 0 the optimization converged for that conformer.
)DOC");

  m.def("UFFGetMoleculeForceField", RDKit::UFFGetMoleculeForceField, "mol"_a,
        "vdwThresh"_a = 10.0, "confId"_a = -1,
        "ignoreInterfragInteractions"_a = true,
        R"DOC(returns a UFF force field for a molecule

ARGUMENTS:

 - mol : the molecule of interest
 - vdwThresh : used to exclude long-range van der Waals interactions
               (defaults to 10.0)
 - confId : indicates which conformer to optimize
 - ignoreInterfragInteractions : if true, nonbonded terms between
               fragments will not be added to the forcefield.
)DOC");

  m.def(
      "UFFHasAllMoleculeParams", RDKit::UFFHasAllMoleculeParams, "mol"_a,
      R"DOC(checks if UFF parameters are available for all of a molecule's atoms

ARGUMENTS:

 - mol : the molecule of interest.
)DOC");

  m.def("MMFFOptimizeMolecule", RDKit::MMFFOptimizeMolecule, "mol"_a,
        "mmffVariant"_a = "MMFF94", "maxIters"_a = 200,
        "nonBondedThresh"_a = 100.0, "confId"_a = -1,
        "ignoreInterfragInteractions"_a = true,
        R"DOC(uses MMFF to optimize a molecule's structure

ARGUMENTS:

 - mol : the molecule of interest
 - mmffVariant : "MMFF94" or "MMFF94s"
 - maxIters : the maximum number of iterations (defaults to 200)
 - nonBondedThresh : used to exclude long-range non-bonded
                interactions (defaults to 100.0)
 - confId : indicates which conformer to optimize
 - ignoreInterfragInteractions : if true, nonbonded terms between
                fragments will not be added to the forcefield

RETURNS: 0 if the optimization converged, -1 if the forcefield could
         not be set up, 1 if more iterations are required.
)DOC");

  m.def("MMFFSanitizeMolecule", RDKit::SanitizeMMFFMol, "mol"_a,
        R"DOC(sanitizes a molecule according to MMFF requirements.

 - mol : the molecule of interest.
)DOC");

  m.def("MMFFGetMoleculeProperties", RDKit::GetMMFFMolProperties, "mol"_a,
        "mmffVariant"_a = "MMFF94", "mmffVerbosity"_a = 0,
        R"DOC(returns an MMFFMolProperties object for a
molecule, which is required by MMFFGetMoleculeForceField()
and can be used to get/set MMFF properties

ARGUMENTS:

 - mol : the molecule of interest
 - mmffVariant : "MMFF94" or "MMFF94s"
               (defaults to "MMFF94")
 - mmffVerbosity : 0: none; 1: low; 2: high (defaults to 0).
)DOC");
  m.def("MMFFGetMoleculeForceField", RDKit::MMFFGetMoleculeForceField, "mol"_a,
        "MMFFMolProperties"_a.none(), "nonBondedThresh"_a = 100.0,
        "confId"_a = -1, "ignoreInterfragInteractions"_a = true,
        R"DOC(returns a MMFF force field for a molecule

ARGUMENTS:

 - mol : the molecule of interest
 - MMFFMolProperties : MMFFMolProperties object as returned
               by MMFFGetMoleculeProperties()
 - nonBondedThresh : used to exclude long-range non-bonded
               interactions (defaults to 100.0)
 - confId : indicates which conformer to optimize
 - ignoreInterfragInteractions : if true, nonbonded terms between
               fragments will not be added to the forcefield
)DOC");

  m.def(
      "CreateEmptyForceFieldForMol", RDKit::CreateEmptyForceFieldForMol,
      "mol"_a, "confId"_a = -1,
      R"DOC(Get An empty Force Field, with only the positions of the atoms but no Contributions.

ARGUMENTS :

    - mol : the molecule of interest
    - confId: the conformer which positions should be added to the force field.
)DOC");

  m.def(
      "MMFFHasAllMoleculeParams", RDKit::MMFFHasAllMoleculeParams, "mol"_a,
      R"DOC(checks if MMFF parameters are available for all of a molecule's atoms

ARGUMENTS:

 - mol : the molecule of interest
)DOC");

  m.def("MMFFOptimizeMoleculeConfs", RDKit::MMFFConfsHelper, "mol"_a,
        "numThreads"_a = 1, "maxIters"_a = 200, "mmffVariant"_a = "MMFF94",
        "nonBondedThresh"_a = 100.0, "ignoreInterfragInteractions"_a = true,
        R"DOC(uses MMFF to optimize all of a molecule's conformations

ARGUMENTS:

 - mol : the molecule of interest
 - numThreads : the number of threads to use, only has an effect if the RDKit
                was built with thread support (defaults to 1)
                If set to zero, the max supported by the system will be used.
 - maxIters : the maximum number of iterations (defaults to 200)
 - mmffVariant : "MMFF94" or "MMFF94s"
 - nonBondedThresh : used to exclude long-range non-bonded
               interactions (defaults to 100.0)
 - ignoreInterfragInteractions : if true, nonbonded terms between
               fragments will not be added to the forcefield.

RETURNS: a list of (not_converged, energy) 2-tuples.
    If not_converged is 0 the optimization converged for that conformer.
)DOC");

  m.def(
      "GetUFFBondStretchParams", ForceFields::getUFFBondStretchParams, "mol"_a,
      "idx1"_a, "idx2"_a,
      "Retrieves UFF bond stretch parameters for atoms with indexes idx1, idx2 "
      "as a (kb, r0) tuple, or None if no parameters could be found");

  m.def("GetUFFAngleBendParams", ForceFields::getUFFAngleBendParams, "mol"_a,
        "idx1"_a, "idx2"_a, "idx3"_a,
        "Retrieves UFF angle bend parameters for atoms with indexes "
        "idx1, idx2, idx3 as a (ka, theta0) tuple, or None if no "
        "parameters could be found");

  m.def("GetUFFTorsionParams", ForceFields::getUFFTorsionParams, "mol"_a,
        "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a,
        "Retrieves UFF torsion parameters for atoms "
        "with indexes idx1, idx2, idx3, idx4 as a V float value, or None "
        "if no parameters could be found");

  m.def("GetUFFInversionParams", ForceFields::getUFFInversionParams, "mol"_a,
        "idx1"_a, "idx2"_a, "idx3"_a, "idx4"_a,
        "Retrieves UFF inversion parameters for atoms "
        "with indexes idx1, idx2, idx3, idx4 as a K float value, or None "
        "if no parameters could be found");

  m.def(
      "GetUFFVdWParams", ForceFields::getUFFVdWParams, "mol"_a, "idx1"_a,
      "idx2"_a,
      "Retrieves UFF van der Waals parameters for atoms with indexes idx1, "
      "idx2 as a (x_ij, D_ij) tuple, or None if no parameters could be found");

  m.def("OptimizeMolecule", RDKit::FFHelper, "ff"_a, "maxIters"_a = 200,
        R"DOC(uses the supplied force field to optimize a molecule's structure

ARGUMENTS:

 - ff : the force field
 - maxIters : the maximum number of iterations (defaults to 200)

RETURNS: 0 if the optimization converged, 1 if more iterations are required.
)DOC");

  m.def(
      "OptimizeMoleculeConfs", RDKit::FFConfsHelper, "mol"_a, "ff"_a,
      "numThreads"_a = 1, "maxIters"_a = 200,
      R"DOC(uses the supplied force field to optimize all of a molecule's conformations

ARGUMENTS:

 - mol : the molecule of interest
 - ff : the force field
 - numThreads : the number of threads to use, only has an effect if the RDKit
                was built with thread support (defaults to 1)
                If set to zero, the max supported by the system will be used.
 - maxIters : the maximum number of iterations (defaults to 200)

RETURNS: a list of (not_converged, energy) 2-tuples.
    If not_converged is 0 the optimization converged for that conformer.
)DOC");
}
