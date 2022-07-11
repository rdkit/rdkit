//
// Copyright (C) 2018-2020 Greg Landrum
//
#include "EHTTools.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif
#include <fstream>

extern "C" {
#include <yaehmop/tightbind/bind.h>
}

namespace RDKit {
namespace EHTTools {
const std::string _EHTCharge = "_EHTCharge";
const std::string _EHTMullikenOP = "_EHTMullikenOP";
const std::string _EHTChargeMatrix = "_EHTChargeMatrix";

// we should only call into the C code, which uses tons of globals, from one
// thread at a time. This mutex enforces that.
#ifdef RDK_BUILD_THREADSAFE_SSS
std::mutex yaehmop_mutex;
#endif

bool runMol(const ROMol &mol, EHTResults &results, int confId,
            bool preserveHamiltonianAndOverlapMatrices) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(yaehmop_mutex);
#endif

  // -----------------------------
  // ----- BOILERPLATE -----------
  // -----------------------------
  FILE *nullfile = fopen("nul", "w");
  status_file = nullfile;
  output_file = nullfile;
  FILE *dest = fopen("run.out", "w+");

  unit_cell = (cell_type *)calloc(1, sizeof(cell_type));
  details = (detail_type *)calloc(1, sizeof(detail_type));

  set_details_defaults(details);
  set_cell_defaults(unit_cell);

  safe_strcpy(details->title, (char *)"RDKit job");

  // molecular calculation
  details->Execution_Mode = MOLECULAR;
  details->num_KPOINTS = 1;
  details->K_POINTS = (k_point_type *)calloc(1, sizeof(k_point_type));
  details->K_POINTS[0].weight = 1.0;
  details->avg_props = 1;
  details->use_symmetry = 1;
  details->find_princ_axes = 0;
  details->net_chg_PRT = 1;
  details->ROP_mat_PRT = 1;
  details->Rchg_mat_PRT = 1;

  unit_cell->using_Zmat = 0;
  unit_cell->using_xtal_coords = 0;
  // ---------------------------------
  // ----- END BOILERPLATE -----------
  // ---------------------------------

  unit_cell->num_atoms = mol.getNumAtoms();
  unit_cell->num_raw_atoms = unit_cell->num_atoms;

  unit_cell->atoms =
      (atom_type *)calloc(unit_cell->num_atoms, sizeof(atom_type));

  const Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    safe_strcpy(unit_cell->atoms[i].symb,
                (char *)mol.getAtomWithIdx(i)->getSymbol().c_str());
    auto p = conf.getAtomPos(i);
    unit_cell->atoms[i].loc.x = p.x;
    unit_cell->atoms[i].loc.y = p.y;
    unit_cell->atoms[i].loc.z = p.z;
  }

  unit_cell->charge = MolOps::getFormalCharge(mol);

  // -----------------------------
  // ----- BOILERPLATE -----------
  // -----------------------------
  const char *parmFilePtr = nullptr;
  std::string pfName = "";
  if (std::getenv("BIND_PARM_FILE") == nullptr) {
    auto rdbase = std::getenv("RDBASE");
    if (rdbase != nullptr) {
      pfName += rdbase;
      pfName += "/Data/eht_parms.dat";
      std::ifstream f(pfName.c_str());
      if (f.good()) {
        parmFilePtr = pfName.c_str();
      } else {
        std::cerr << "file " << pfName << " doesn't seem to exist" << std::endl;
      }
    }
  }

  fill_atomic_parms(unit_cell->atoms, unit_cell->num_atoms, nullptr,
                    const_cast<char *>(parmFilePtr));
  unit_cell->num_raw_atoms = unit_cell->num_atoms;
  charge_to_num_electrons(unit_cell);
  build_orbital_lookup_table(unit_cell, &num_orbs, &orbital_lookup_table);

  run_eht(dest);
  // ---------------------------------
  // ----- END BOILERPLATE -----------
  // ---------------------------------

  // pull properties
  results.numAtoms = mol.getNumAtoms();
  results.numOrbitals = num_orbs;
  results.numElectrons = std::lround(unit_cell->num_electrons);
  results.fermiEnergy = properties.Fermi_E;
  results.totalEnergy = properties.total_E;
  results.atomicCharges = std::make_unique<double[]>(mol.getNumAtoms());
  std::memcpy(static_cast<void *>(results.atomicCharges.get()),
              static_cast<void *>(properties.net_chgs),
              mol.getNumAtoms() * sizeof(double));
  size_t sz = mol.getNumAtoms() * num_orbs;
  results.reducedChargeMatrix = std::make_unique<double[]>(sz);
  memcpy(static_cast<void *>(results.reducedChargeMatrix.get()),
         static_cast<void *>(properties.Rchg_mat), sz * sizeof(double));

  sz = mol.getNumAtoms() * (mol.getNumAtoms() + 1) / 2;
  results.reducedOverlapPopulationMatrix = std::make_unique<double[]>(sz);
  memcpy(static_cast<void *>(results.reducedOverlapPopulationMatrix.get()),
         static_cast<void *>(properties.ROP_mat), sz * sizeof(double));

  results.orbitalEnergies = std::make_unique<double[]>(num_orbs);
  std::memcpy(static_cast<void *>(results.orbitalEnergies.get()),
              static_cast<void *>(eigenset.val), num_orbs * sizeof(double));

  if (preserveHamiltonianAndOverlapMatrices) {
    // these need to be recalculated, because they were overwritten during the
    // calculation
    R_space_overlap_matrix(unit_cell, details, Overlap_R, num_orbs,
                           tot_overlaps, orbital_lookup_table, 0);
    full_R_space_Hamiltonian(unit_cell, details, Overlap_R, Hamil_R, num_orbs,
                             orbital_lookup_table, 1);
    sz = num_orbs * num_orbs * sizeof(double);
    results.hamiltonianMatrix = std::make_unique<double[]>(sz);
    std::memcpy(static_cast<void *>(results.hamiltonianMatrix.get()),
                static_cast<void *>(Hamil_R.mat), sz);
    results.overlapMatrix = std::make_unique<double[]>(sz);
    std::memcpy(static_cast<void *>(results.overlapMatrix.get()),
                static_cast<void *>(Overlap_R.mat), sz);
  }
  cleanup_memory();

  fclose(nullfile);
  fclose(dest);
  return true;
}

}  // end of namespace EHTTools
}  // end of namespace RDKit