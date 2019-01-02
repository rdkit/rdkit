//
// Copyright (C) 2018 Greg Landrum
//
#include "EHTTools.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <mutex>


extern "C" {
//#include "local.h"
#include <bind.h>
}

namespace RDKit {
namespace EHTTools {
const std::string _EHTCharge = "_EHTCharge";
const std::string _EHTMullikenOP = "_EHTMullikenOP";
const std::string _EHTChargeMatrix = "_EHTChargeMatrix";


// we should only call into the C code, which uses tons of globals, from one thread
// at a time. This mutex enforces that.
std::mutex yaehmop_mutex;

bool runMol(const ROMol &mol, int confId){
  std::lock_guard<std::mutex> lock(yaehmop_mutex);

  // -----------------------------
  // ----- BOILERPLATE -----------
  // -----------------------------
  FILE *nullfile = fopen("nul", "w");
  status_file = nullfile;
  output_file = nullfile;
  FILE *dest=fopen("run.out","w+");

  unit_cell = (cell_type *)calloc(1,sizeof(cell_type));
  details = (detail_type *)calloc(1, sizeof(detail_type));

  /*******
    initialize some variables
    ********/
  details->walsh_details.num_steps = 1;
  details->walsh_details.num_vars = 0;
  details->use_symmetry = 0;
  details->find_princ_axes = 0;
  details->vary_zeta = 0;
  details->avg_props = 0;
  details->just_geom = 0;
  details->dump_overlap = 0;
  details->dump_hamil = 0;
  details->sparsify_value = 0.0;
  details->Execution_Mode = FAT;
  details->the_const = THE_CONST;
  details->weighted_Hij = 1;
  details->eval_electrostat = 0;
  details->close_nn_contact = NN_DISTANCE;
  details->symm_tol = 0.01;
  details->muller_mix = MULLER_MIX_DEF;
  details->muller_E_tol = MULLER_E_TOL_DEF;
  details->muller_Z_tol = MULLER_Z_TOL_DEF;
  details->num_moments = 4;
  details->line_width = 80;
  details->k_offset = K_OFFSET;
  unit_cell->equiv_atoms = 0;


  safe_strcpy(details->title,(char *)"RDKit job");


  // molecular calculation
  details->Execution_Mode = MOLECULAR;
  details->num_KPOINTS = 1;
  details->K_POINTS = (k_point_type *)calloc(1, sizeof(k_point_type));
  details->K_POINTS[0].weight = 1.0;
  details->avg_props = 0;
  details->use_symmetry = 1;
  details->find_princ_axes = 0;
  details->net_chg_PRT=1;
  details->ROP_mat_PRT=1;
  details->Rchg_mat_PRT=1;

  unit_cell->using_Zmat = 0;
  unit_cell->using_xtal_coords = 0;
  // ---------------------------------
  // ----- END BOILERPLATE -----------
  // ---------------------------------

  unit_cell->num_atoms = mol.getNumAtoms();
  unit_cell->num_raw_atoms = unit_cell->num_atoms;

  unit_cell->atoms = (atom_type *)calloc(unit_cell->num_atoms, sizeof(atom_type));

  const Conformer &conf = mol.getConformer(confId);
  for(unsigned int i=0;i<mol.getNumAtoms();++i){
    safe_strcpy(unit_cell->atoms[i].symb,(char *)mol.getAtomWithIdx(i)->getSymbol().c_str());
    auto p = conf.getAtomPos(i);
    unit_cell->atoms[i].loc.x = p.x;
    unit_cell->atoms[i].loc.y = p.y;
    unit_cell->atoms[i].loc.z = p.z;
  }

  unit_cell->charge = MolOps::getFormalCharge(mol);

  // -----------------------------
  // ----- BOILERPLATE -----------
  // -----------------------------
  fill_atomic_parms(unit_cell->atoms,unit_cell->num_atoms,NULL,NULL);
  unit_cell->num_raw_atoms = unit_cell->num_atoms;
  charge_to_num_electrons(unit_cell);
  build_orbital_lookup_table(unit_cell,&num_orbs,&orbital_lookup_table);

  run_eht(dest);
  // ---------------------------------
  // ----- END BOILERPLATE -----------
  // ---------------------------------

  //pull properties
  for(auto &atom : mol.atoms() ){
    atom->setProp(_EHTCharge,properties.net_chgs[atom->getIdx()],true);
  }

  size_t sz = mol.getNumAtoms()*mol.getNumAtoms();
  auto ropMat = new double[sz];
  memcpy(ropMat,properties.ROP_mat,sz*sizeof(double));
  boost::shared_array<double> sptr(ropMat);
  mol.setProp(_EHTMullikenOP,sptr,true);
  auto rchgMat = new double[sz];
  memcpy(rchgMat,properties.Rchg_mat,sz*sizeof(double));
  boost::shared_array<double> sptr2(rchgMat);
  mol.setProp(_EHTChargeMatrix,sptr2,true);

  cleanup_memory();

  fclose(nullfile);
  fclose(dest);
  return true;
}


} // end of namespace EHTTools
} // end of namespace RDKit