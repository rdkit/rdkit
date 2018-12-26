//
// Copyright (C) 2018 Greg Landrum
//
#include "EHTTools.h"
#include <GraphMol/RDKitBase.h>

extern "C" {
//#include "local.h"
#include <bind.h>
}


namespace EHTTools {
  using namespace RDKit;


void stub() {
  char err_string[240];

  FILE *nullfile = fopen("nul","w");

  status_file = nullfile;
  output_file = nullfile;

  unit_cell = (cell_type *)calloc(1,sizeof(cell_type));
  details = (detail_type *)calloc(1,sizeof(detail_type));
  if(!unit_cell || !details) fatal("Can't allocate initial memory.");
  FILE *dest=nullfile;

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
  details->symm_tol = SYMM_TOL;
  details->muller_mix = MULLER_MIX_DEF;
  details->muller_E_tol = MULLER_E_TOL_DEF;
  details->muller_Z_tol = MULLER_Z_TOL_DEF;
  details->num_moments = 4;
  details->line_width = 80;
  details->k_offset = K_OFFSET;
  unit_cell->equiv_atoms = 0;


  safe_strcpy(details->title,"test job");


  // molecular calculation
  details->Execution_Mode = MOLECULAR;
  details->num_KPOINTS = 1;
  details->K_POINTS = (k_point_type *)calloc(1,sizeof(k_point_type));
  if( !details->K_POINTS ) fatal("Can't allocate the single k point.");
  details->K_POINTS[0].weight = 1.0;
  details->avg_props = 0;
  details->use_symmetry = 1;
  details->net_chg_PRT=1;
  details->ROP_mat_PRT=1;


  unit_cell->using_Zmat = 0;
  unit_cell->using_xtal_coords = 0;


  unit_cell->num_atoms = 3;
  unit_cell->atoms = (atom_type *)calloc(unit_cell->num_atoms,sizeof(atom_type));
  if(!unit_cell->atoms){
    sprintf(err_string,"Can't allocate memory for: %d atoms.",unit_cell->num_atoms);
    fatal("Can't allocate memory for the atoms.");
  }
  safe_strcpy(unit_cell->atoms[0].symb,"H");
  unit_cell->atoms[0].loc.x = 0.0;
  unit_cell->atoms[0].loc.y = 0.0;
  unit_cell->atoms[0].loc.z = 0.0;

  safe_strcpy(unit_cell->atoms[1].symb,"C");
  unit_cell->atoms[1].loc.x = 0.0;
  unit_cell->atoms[1].loc.y = 0.0;
  unit_cell->atoms[1].loc.z = 1.0;

  safe_strcpy(unit_cell->atoms[2].symb,"N");
  unit_cell->atoms[2].loc.x = 0.0;
  unit_cell->atoms[2].loc.y = 0.0;
  unit_cell->atoms[2].loc.z = 2.1;

  unit_cell->charge = 0.0;


  // shouldn't need to change below here
  fill_atomic_parms(unit_cell->atoms,unit_cell->num_atoms,NULL,NULL);
  unit_cell->num_raw_atoms = unit_cell->num_atoms;
  charge_to_num_electrons(unit_cell);
  build_orbital_lookup_table(unit_cell,&num_orbs,&orbital_lookup_table);


  /* install the sig_int handler */
  //signal(SIGINT,handle_sigint);

  run_eht(dest);

  //pull properties
  for(int i=0;i<unit_cell->num_atoms;i++){
   printf(">>>> Atom %d: %.2f\n",i+1,properties.net_chgs[i]);
  }

  for(int i=0;i<unit_cell->num_atoms;i++){
  for(int j=0;j<i;j++){
   printf(">>>> ROP %d-%d: %.2f\n",i+1,j+1,properties.ROP_mat[i*(i+1)/2 + j]);
  }
}

  free(unit_cell);
  free(details);

}



}
