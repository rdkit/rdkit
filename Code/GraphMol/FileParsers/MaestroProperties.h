//
//  Copyright (C) 2024 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>

namespace RDKit {
namespace FileParsers {
namespace schrodinger {

static const std::string MAE_ATOM_TYPE = "i_m_mmod_type";
static const std::string MAE_BOND_DATIVE_MARK = "b_sPrivate_dative_bond";
static const std::string MAE_BOND_PARITY = "i_sd_original_parity";
static const std::string MAE_ENHANCED_STEREO_STATUS =
    "i_m_ct_enhanced_stereo_status";
static const std::string MAE_RGROUP_LABEL = "i_sd__MolFileRLabel";
static const std::string MAE_STEREO_STATUS = "i_m_ct_stereo_status";
static const std::string PDB_ATOM_NAME = "s_m_pdb_atom_name";
static const std::string PDB_CHAIN_NAME = "s_m_chain_name";
static const std::string PDB_INSERTION_CODE = "s_m_insertion_code";
static const std::string PDB_OCCUPANCY = "r_m_pdb_occupancy";
static const std::string PDB_RESIDUE_NAME = "s_m_pdb_residue_name";
static const std::string PDB_RESIDUE_NUMBER = "i_m_residue_number";
static const std::string PDB_TFACTOR = "r_m_pdb_tfactor";

}  // namespace schrodinger
}  // namespace FileParsers
}  // namespace RDKit
