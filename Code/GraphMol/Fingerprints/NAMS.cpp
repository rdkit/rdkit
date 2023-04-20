//
//  Copyright (C) 2023 Rocco Moretti and other RDKit contributors
//  Copyright (C) 2013, Andre Falcao and Ana Teixeira, University of Lisbon - LaSIGE
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Portions of this implementation include or are based on NAMS by Andre Falcao and Ana Teixeira
//  (https://github.com/aofalcao/nams-docker)
//  NAMS is free software: it can be redistributed and/or modifed
//  under the terms of the MIT License as published on the official site of Open Source Initiative
//
// Please cite the authors in any work or product based on this material:
//
// AL Teixeira, AO Falcao. 2013. A non-contiguous atom matching structural similarity function. J. Chem. Inf. Model. DOI: 10.1021/ci400324u.
// (http://dx.doi.org/10.1021/ci400324u)

#include <GraphMol/Fingerprints/NAMS.h>
#include <GraphMol/Fingerprints/hungarian.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>

#include <GraphMol/SmilesParse/SmilesWrite.h>
//#include <GraphMol/Descriptors/MolDescriptors.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <tuple>
#include <map>

namespace RDKit {
namespace NAMS {

const int LU_ELEMS_DEFAULT[] ={0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 8, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static constexpr int NUM_ELEMS = 11;

const float ADM0[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f  } };

const float ADM1[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f  } };

const float ADM2[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  0.48f ,  0.64f ,  0.80f ,  0.96f ,  0.65f ,  0.81f ,  0.97f ,  0.98f ,  0.98f ,  1.00f  },
 { 0.48f,  0.00f ,  0.16f ,  0.50f ,  0.48f ,  0.70f ,  0.50f ,  0.48f ,  0.90f ,  0.50f ,  0.52f  },
 { 0.64f,  0.16f ,  0.00f ,  0.20f ,  0.50f ,  0.40f ,  0.80f ,  0.60f ,  0.70f ,  0.70f ,  0.80f  },
 { 0.80f,  0.50f ,  0.20f ,  0.00f ,  0.16f ,  0.17f ,  0.90f ,  0.17f ,  0.01f ,  0.21f ,  0.27f  },
 { 0.96f,  0.48f ,  0.50f ,  0.16f ,  0.00f ,  0.33f ,  0.17f ,  0.07f ,  0.90f ,  0.14f ,  0.21f  },
 { 0.65f,  0.70f ,  0.40f ,  0.17f ,  0.33f ,  0.00f ,  0.50f ,  0.32f ,  0.07f ,  0.33f ,  0.35f  },
 { 0.81f,  0.50f ,  0.80f ,  0.90f ,  0.17f ,  0.50f ,  0.00f ,  0.16f ,  0.80f ,  0.17f ,  0.21f  },
 { 0.97f,  0.48f ,  0.60f ,  0.17f ,  0.07f ,  0.32f ,  0.16f ,  0.00f ,  0.85f ,  0.07f ,  0.14f  },
 { 0.98f,  0.90f ,  0.70f ,  0.01f ,  0.90f ,  0.07f ,  0.80f ,  0.85f ,  0.00f ,  0.80f ,  0.90f  },
 { 0.98f,  0.50f ,  0.70f ,  0.21f ,  0.14f ,  0.33f ,  0.17f ,  0.07f ,  0.80f ,  0.00f ,  0.07f  },
 { 1.00f,  0.52f ,  0.80f ,  0.27f ,  0.21f ,  0.35f ,  0.21f ,  0.14f ,  0.90f ,  0.07f ,  0.00f  } };

const float ADM3[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.00f ,  0.80f ,  0.90f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.80f ,  0.00f ,  1.00f ,  1.00f ,  0.70f ,  1.00f ,  1.00f ,  0.90f ,  1.00f ,  1.00f  },
 { 1.00f,  0.90f ,  1.00f ,  0.00f ,  0.80f ,  1.00f ,  0.50f ,  0.90f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.80f ,  0.00f ,  1.00f ,  1.00f ,  0.10f ,  1.00f ,  0.20f ,  0.30f  },
 { 1.00f,  1.00f ,  0.70f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.50f ,  1.00f ,  1.00f ,  0.00f ,  0.90f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.90f ,  0.10f ,  1.00f ,  0.90f ,  0.00f ,  1.00f ,  0.10f ,  0.20f  },
 { 1.00f,  1.00f ,  0.90f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.20f ,  1.00f ,  1.00f ,  0.10f ,  1.00f ,  0.00f ,  0.10f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.30f ,  1.00f ,  1.00f ,  0.20f ,  1.00f ,  0.10f ,  0.00f  } };

const float ADM4[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  } };


NAMSParameters::~NAMSParameters() {
  free(blev_mat);
}

float
NAMSParameters::ELEMS_DISTS(int ele1, int ele2) const {
  int a1 = LU_ELEMS_DEFAULT[ele1];
  int a2 = LU_ELEMS_DEFAULT[ele2];

  switch( ADM ) {
  case 0:
    return ADM0[a1][a2];
  case 1:
    return ADM1[a1][a2];
  case 2:
    return ADM2[a1][a2];
  case 3:
    return ADM3[a1][a2];
  default:
    return ADM4[a1][a2];
  }
}

const int*
NAMSParameters::getBondLevelsMatrix() const {
  if ( blev_mat == nullptr || blev_alpha != BS_ALPHA ) {
    calcBondLevelsMatrix();
  }
  return blev_mat;
}

// From nams-docker/nams/nams.cpp:calcBondLevelsMatrix()
void
NAMSParameters::calcBondLevelsMatrix() const {

  //this is a very important function and for efficiency should be computed only ONCE and the results stored in one array
  //the idea is that we should compute the floating point stuf once and only once
  // to make everything really fast, the striangular matrix will be symmetrical
  float alpha = this->BS_ALPHA;
  int *mat=(int *)malloc(sizeof(int)*MAX_LEVELS*MAX_LEVELS);
  int v, mx;
  for(int i=0; i<MAX_LEVELS-1; i++) {
    for(int j=i; j<MAX_LEVELS; j++) {
      //v=(int)(100*1.0f/powf((abs(i-j)+j+1.0f), parms->BS_ALPHA));
      //modified at 20/11/2017
      mx=i;
      if(alpha>0.0f) {
        if(j>mx) mx=j;
        v=(int)(100.0*powf(alpha, (float)(mx+abs(i-j))));
      }else {
        if(i==j) v=100; else v=0;
      }
      mat[i*MAX_LEVELS+j] = v;
      mat[j*MAX_LEVELS+i] = v;
    }
  }

  blev_mat = mat;
  blev_alpha = alpha;
}

/*!
  From nams-docker/makenamsdb/recoder.py:Recoder::process_bonds()

  #the idea is to select all the immediate bonds in the start set
  # append them to the current level_bonds and elliminate them from the connection matrix

  Recursive algorithm.

  \param mol: The molecule being processed
  \param start_set: The atom indexes to process for this level
  \param eliminated: A set of atom pairs (bonds) which have already been processed (implementation detail
  \param lvl_bonds: A vector of levels, where each level contains a vector of atom pairings (return-by-reference)

*/
void process_bonds(
  const ROMol &mol,
  const std::vector< unsigned int > & start_set,
  std::set< std::pair< unsigned int, unsigned int > > & eliminated,
  std::vector< std::vector< std::pair< unsigned int, unsigned int > > > & lvl_bonds
) {

  std::vector< unsigned int > to_follow;
  std::vector< std::pair< unsigned int, unsigned int > > this_level;
  for ( unsigned int at1: start_set ) {
    for ( const Atom* atom2: mol.atomNeighbors(mol.getAtomWithIdx(at1)) ) {
      unsigned int at2 = atom2->getIdx();
      std::pair< unsigned int, unsigned int > this_bond{ at1, at2 };
      if ( eliminated.count( this_bond ) > 0 ) { continue; }

      this_level.push_back( this_bond );
      eliminated.insert( this_bond );
      eliminated.insert( std::make_pair( at2, at1 ) );
      to_follow.push_back( at2 );
    }
  }
  lvl_bonds.push_back( this_level );

  if ( to_follow.size() > 0 ) {
    process_bonds( mol, to_follow, eliminated, lvl_bonds );
  }
}

BondInfoType::BondInfoType(const ROMol &mol, unsigned int at1, unsigned int at2, bool do_isomerism) {
  const Atom* atom1 = mol.getAtomWithIdx(at1);
  const Atom* atom2 = mol.getAtomWithIdx(at2);
  ele1 = atom1->getAtomicNum();
  ele2 = atom2->getAtomicNum();
  nrings1 = mol.getRingInfo()->numAtomRings(at1);
  nrings2 = mol.getRingInfo()->numAtomRings(at2);
  if ( do_isomerism ) { /*TODO: implement chirality?*/ }
  const Bond* bond = mol.getBondBetweenAtoms(at1,at2);
  assert( bond != nullptr );
  inring = mol.getRingInfo()->minBondRingSize(bond->getIdx()) != 0;
  aromatic = bond->getIsAromatic();
  order = bond->getBondTypeAsDouble();
  if ( do_isomerism ) { /*TODO: implement isomerism?*/ }
}

bool
BondInfoType::operator<( const BondInfoType & r ) const {
  return std::make_tuple(ele1, nrings1, chir1, ele2, nrings2, chir2, inring, aromatic, order, dbcistrans)
    < std::make_tuple(r.ele1, r.nrings1, r.chir1, r.ele2, r.nrings2, r.chir2, r.inring, r.aromatic, r.order, r.dbcistrans);

}

class BondInfo {
public:
  BondInfo(const ROMol &mol, unsigned int at1, unsigned int at2, bool do_isomerism = false):
    atm1(at1),
    atm2(at2),
    type(mol, at1, at2, do_isomerism)
  {}

  unsigned int atm1, atm2;
  BondInfoType type;
};

/*!
 From nams-docker/makenamsdb/recoder.py:Recoder::get_mol_info()

 Returns a pre-processed data structure for the molecule
 Top level is a vector with one entry per atom
 Second level is a vector with one entry per distance step
 Third level is a vector of bond information for all the bonds at that level
*/
std::vector< std::vector< std::vector< BondInfo > > >
get_connection_matrix(const ROMol &mol, bool do_isomerism = false) {

  MolOps::findSSSR(mol); // Add ring information to mol's RingInfo object

  std::vector< std::vector< std::vector< BondInfo > > > connection_matrix;

  for( unsigned int atomno = 0; atomno < mol.getNumAtoms(); ++atomno ) {
    std::vector< unsigned int > start_set;
    start_set.push_back( atomno );

    std::set< std::pair< unsigned int, unsigned int > > eliminated;
    std::vector< std::vector< std::pair< unsigned int, unsigned int > > > lvl_bonds;
    process_bonds( mol, start_set, eliminated, lvl_bonds );

    std::vector< std::vector< BondInfo > > bonds_for_all_levels;
    for ( auto const & level: lvl_bonds ) {
      if ( level.empty() ) { continue; } // "for some weird reason the last level is always empty so take it out"
      std::vector< BondInfo > bonds_for_level;
      for ( auto const & bond_pair: level ) {
        bonds_for_level.emplace_back( mol, bond_pair.first, bond_pair.second, do_isomerism );
      }
      bonds_for_all_levels.emplace_back( std::move(bonds_for_level) );
    }

    connection_matrix.emplace_back( std::move(bonds_for_all_levels) );
   }

  return connection_matrix;
}

NAMSMolInfo::NAMSMolInfo(const ROMol &mol) {
  smiles = MolToSmiles( mol );
//  molwt = Descriptors::calcAMW(mol);

  std::vector< std::vector< std::vector< BondInfo > > > connection_matrix = get_connection_matrix(mol);
  // From nams-docker/makenamsdb/recoder.py:Recoder::export_mol_info()

  unsigned int natoms = connection_matrix.size();
  mat_aba_types.resize( natoms );
  mat_levels.resize( natoms );
  // Mapping from the BondInfoType to the index
  std::map< BondInfoType, unsigned int > aba_type_index;

  for ( unsigned int atmno=0; atmno < natoms; ++atmno ) {
    auto const & atominfo = connection_matrix[atmno];
    unsigned int nlevels = atominfo.size();
    for( unsigned int L=0; L < nlevels; ++L ) {
      for ( const BondInfo & binfo: atominfo[L] ) {
        auto itr = aba_type_index.find(binfo.type);
        if ( itr == aba_type_index.end() ) {
          aba_type_index[ binfo.type ] = aba_type_index.size(); // 0 indexed
          itr = aba_type_index.find(binfo.type);
        }
        mat_aba_types[atmno].push_back( itr->second );
        mat_levels[atmno].push_back( L );
      }
    }
  }

  aba_types.resize( aba_type_index.size() );
  for ( auto const & entry: aba_type_index ) {
    aba_types[ entry.second ] = entry.first;
  }

}

unsigned int
NAMSMolInfo::natoms() const {
  assert( mat_aba_types.size() == mat_levels.size() );
  return mat_aba_types.size();
}

unsigned int
NAMSMolInfo::nbonds() const {
  if ( mat_aba_types.empty() ) { return 0; }
  return mat_aba_types[0].size();
}

unsigned int
NAMSMolInfo::naba_types() const {
  return aba_types.size();
}

std::ostream & operator<<( std::ostream & os, BondInfoType const & type ) {
  os << type.ele1 << ' ' << type.nrings1 << ' ' << type.chir1;
  os << ' ' << type.ele2 << ' ' << type.nrings2 << ' ' << type.chir2;
  os << ' ' << type.inring << ' ' << type.aromatic << ' ' << type.order*10 << ' ' << type.dbcistrans;
  os << '\n';
  return os;
}

// This should dump information to cout in approx
std::string
NAMSMolInfo::dump(int cid) const {
  std::stringstream ss;

  unsigned int nbonds = 0;
  if ( ! mat_aba_types.empty() ) {
    nbonds = mat_aba_types[0].size();
  }

  ss << cid << ' ' << int(molwt*10) << ' ' << smiles << '\n';
  ss << natoms() << ' ' << nbonds << ' ' << aba_types.size() << '\n';
  for ( auto const & type: aba_types ) {
    ss << type;
  }
  for ( auto const & bond: mat_aba_types ) {
    for ( unsigned int x: bond ) {
      ss << x << ' ';
    }
    ss << '\n';
  }
  for ( auto const & bond: mat_levels ) {
    for ( unsigned int x: bond ) {
      ss << x << ' ';
    }
    ss << '\n';
  }

  return ss.str();
}

NAMSMolInfo * getNAMSMolInfo(const ROMol &mol) {
  return new NAMSMolInfo(mol);
}

// from nams-docker/nams/nams.cpp:compare_aba_bonds(), with minimal tweaks
int compare_aba_bonds(const BondInfoType & aba1, const BondInfoType & aba2, const NAMSParameters & parms)
{

  float simil=1.0;
  //now we can compute the atom similarity
  simil *= (1.0f - parms.ELEMS_DISTS(aba1.ele1,aba2.ele1));
  simil *= (1.0f - parms.ELEMS_DISTS(aba1.ele2,aba2.ele2));

  if(aba1.chir1*aba2.chir1 ==-1) simil *= parms.ACHIR_FAC;
  if(aba1.chir2*aba2.chir2 ==-1) simil *= parms.ACHIR_FAC;

  simil *= powf(parms.ANRINGS_FAC, (float)std::abs(aba1.nrings1-aba2.nrings1));
  simil *= powf(parms.ANRINGS_FAC, (float)std::abs(aba1.nrings2-aba2.nrings2));

  if(aba1.inring != aba2.inring) simil *= parms.BRING_FAC;
  if(aba1.aromatic != aba2.aromatic) simil   *= parms.BAROM_FAC;
  if(aba1.order != aba2.order) simil  *= parms.BORDER_FAC;
  if(aba1.dbcistrans * aba2.dbcistrans == -1) simil *= parms.DBSTEREO_FAC;

  return (int)(simil*100);
}


// from nams-docker/nams/nams.cpp:getAbaBondsCompMatrix(), with minimal tweaks
// The caller is responsible for calling malloc's free() on the returned value
int *getAbaBondsCompMatrix(const NAMSMolInfo & mi1, const NAMSMolInfo & mi2, const NAMSParameters & parms)
{
  //this function receives two molecules and will produce a matrix with the comparison of all bonds one agains the other
  //using the order of the ababonds in each atom
  // the matrix is represented as an array of integers(?) to access each the formula is
  //mat[a2+na2*a1]
  int *mat=NULL;
  int res;

  mat=(int *)malloc(mi1.naba_types()*mi2.naba_types()*sizeof(int));

  for(unsigned int i=0; i<mi1.naba_types(); i++) {
    for(unsigned int j=0; j<mi2.naba_types(); j++) {
      res=compare_aba_bonds(mi1.aba_types[i], mi2.aba_types[j], parms);
      mat[mi2.naba_types()*i+ j]=res;
    }
  }

  return mat;
}

// from nams-docker/nams/nams.cpp:matchBonds(), with minimal tweaks
int matchBonds(int a1, int a2, const NAMSMolInfo & mi1, const NAMSMolInfo & mi2, int *aba_ms, const NAMSParameters & parms)
{
  // this is one of the MOST CRITICAL functions of this whole story. Must be carefully crafted
  // here we will make a matrix with as many rows as the number of bonds of mol1 and rows as n of bonds for mol2
  // the matrix will have the comparison score USING the levels for each bond of mol1 against each bond of mol2
  // aba_ms -> similarity matrix between aba_bond types
  // blev_mat -> levels matrix - >precomputed factors that account for the level distances between bonds
  // a1 and a2 -> atom indexes for getting the corresponding rows out of the molInfos

  int lev1, lev2, abat1, abat2, nfac = mi2.naba_types();
  int sim=0, mscore=0, ce1=0, ce2=0;
  //thes variables are just for information purposes of counting the best possible matches on rows and columns
  int max_sim=0;
  int max_row=0, max_col=0;

  const int* blev_mat = parms.getBondLevelsMatrix();

  //SimBond *sbs;  //this probably could be allocated globally as for the levels matrix and used over and over again for all molecules
  //bool *elims1, *elims2; //for these two the same caveat applies;
  int nbonds1=mi1.nbonds();
  int nbonds2=mi2.nbonds();

  std::vector< bool > elims1( nbonds1, false );
  std::vector< bool > elims2( nbonds2, false );
  std::vector< std::vector< int > > the_mat;
  the_mat.resize(nbonds1);
  for ( int ii=0; ii<nbonds1; ++ii) {
    the_mat[ii].resize(nbonds2);
  }

  std::vector< unsigned int > const & abas1 = mi1.mat_aba_types[a1];
  std::vector< unsigned int > const & abas2 = mi2.mat_aba_types[a2];
  std::vector< int > const & levs1 = mi1.mat_levels[a1];
  std::vector< int > const & levs2 = mi2.mat_levels[a2];


  //this is it! Compute the bond matrix against 2 atoms
  //printf("bond matrix against 2 atoms %d %d\n",mi1->nbonds, mi2->nbonds );
  for(int i=0; i<nbonds1; i++) {
    abat1 = abas1[i];
    lev1  = levs1[i];
    //printf("----->");
    for(int j=0; j<nbonds2; j++) {
      abat2 = abas2[j];
      lev2  = levs2[j];
      sim=blev_mat[lev1*NAMSParameters::MAX_LEVELS+lev2] * aba_ms[abat1*nfac+abat2];
      if(sim>max_sim) {
        max_sim = sim;
        max_row=i;
        max_col=j;
      }
      //the_mat[i*nbonds2+j] = sim;
      the_mat[i][j]=sim;
    }
  }




  elims1[max_row]=true;
  elims2[max_col]=true;
  mscore=max_sim;
  //memset(the_mat[max_row], -1, sizeof(int)*nbonds2);
  //for(int i=0; i<nbonds1; i++) the_mat[i][max_col]=-1;

  //this could speed up by counting the number of hits and breaking when hits = nbonds
  while(ce1<nbonds1 && ce2 <nbonds2) {
    max_sim=-1;
    for(int row=0; row< nbonds1; row++) {
      if(!elims1[row]) {
        for(int col=0; col<nbonds2; col++) {
          if(!elims2[col] && the_mat[row][col] > max_sim) {
            //max_sim=the_mat[row][col];
            max_sim=the_mat[row][col];
            max_row=row;
            max_col=col;
          }
        }
      }
    }
    elims1[max_row]=true;
    elims2[max_col]=true;
    mscore+=max_sim;
    //memset(the_mat[max_row], -1, sizeof(int)*nbonds2);
    //for(int i=0; i<nbonds1; i++) the_mat[i][max_col]=-1;
    ce1++;
    ce2++;
  }

  return mscore;
}


// from nams-docker/nams/nams.cpp:calcSelfSimilarity(), with minimal tweaks
double calcSelfSimilarity(const NAMSMolInfo & mi, const NAMSParameters & parms)
{
  //this function returns the self similarity in one molecule by simply computing the diagonal elements
  int*  mat=(int*)malloc(mi.natoms()*mi.natoms()*sizeof(int));
  int bscore=0;
  int * aba_match_scores = getAbaBondsCompMatrix(mi, mi, parms);

  for(unsigned int i=0; i<mi.natoms();i++) {
    bscore += matchBonds(i, i, mi, mi, aba_match_scores, parms);
  }
  free(mat);
  free(aba_match_scores);
  return bscore/10000.0d;
}

// from nams-docker/nams/nams.cpp:calcAtomMatchingMatrix(), with minimal tweaks
int *calcAtomMatchingMatrix(const NAMSMolInfo & mi1, const NAMSMolInfo & mi2, const NAMSParameters & parms, int *aba_match_scores)
{
  //this is a critical function. It will compare all atoms on mol1 to all atoms of mol2 and create a score for each match
  //mat will hold the scores
  int*  mat=(int*)malloc(mi1.natoms()*mi2.natoms()*sizeof(int));

  int c=0, bscore;
  for(unsigned int i=0; i<mi1.natoms();i++) {
    //VERY VERY IMPORTANT! WHEN WEIGTHING THE ATOMS IMPORTANCE THIS IS WHERE WE SHOULD ACCOUNT FOR THEM
    for(unsigned int j=0; j < mi2.natoms(); j++) {
      bscore = matchBonds(i, j, mi1, mi2, aba_match_scores, parms);
      //the penalties stuff
      bscore -= (int)(parms.PEN*abs(int(mi1.nbonds())-int(mi2.nbonds()))*100.0);
      mat[c]=bscore;
      c++;
      //printf("%5.2f ", bscore/10000.0f);
    }
    //puts("");
  }
  return mat;
}

// from nams-docker/nams/nams.cpp:nams_runner(), with minimal tweaks
void nams_runner(const NAMSMolInfo & mi1, const NAMSMolInfo & mi2, const NAMSParameters & parms, NAMSResult & result, float *wts)
{
  //Here is the backbone of it all. 2 molecules to be compared
  //d) the atom matching matrix (M)
  //e) the alignment produced (A)

  int  *aba_match_scores; //the scores for each aba_bond match between 2 molecules
  int  *atom_matrix;
  int  v, final_score=0;
  hungarian_problem_t p;

  aba_match_scores = getAbaBondsCompMatrix(mi1, mi2, parms);

  //for(int i=0; i<mi->naba_types; i++) {
  //  for(int j=0; j<mi->naba_types; j++) printf("%6d ", aba_match_scores[mi->naba_types*i+ j]);
  //  puts("");
  //}

  atom_matrix=calcAtomMatchingMatrix(mi1, mi2, parms, aba_match_scores);
  //puts("");


  int** m = array_to_matrix(atom_matrix, mi1.natoms(), mi2.natoms()); //THIS is where the weights will strike in!
  if(wts) {
    for(unsigned int row=0; row<mi1.natoms(); row++) {
      for(unsigned int col=0; col<mi2.natoms(); col++) {
        if(m[row][col]>0) m[row][col]=(int)(powf(m[row][col]/10000.0f, wts[row])*10000);
        else m[row][col] = 0;
      }
    }
  }

  /*int matrix_size =*/ hungarian_init(&p, m , mi1.natoms(), mi2.natoms(), HUNGARIAN_MODE_MAXIMIZE_UTIL) ;

  //here it is!!!
  hungarian_solve(&p);
//  //this is just for writing
//  if(M) {
//    puts("Atom matching matrix");
//    printf("    ");
//    for(unsigned int j =0; j<mi2.natoms(); j++) printf("%6d ", j+1);
//    puts("");
//    for(unsigned int i =0; i<mi1.natoms(); i++) {
//      printf("%3d ", i+1);
//      for(unsigned int j =0; j<mi2.natoms(); j++) {
//        printf("%6.2f ", atom_matrix[i*mi2.natoms()+j]/10000.0f);
//      }
//      puts("");
//    }
//  }
  //this is where the actual result is calculated
  result.mapping1to2.clear();
  result.mapping1to2.resize(mi1.natoms(), -1);
  result.atom_scores.clear();
  result.atom_scores.resize(mi1.natoms(), 0.0);
  for (unsigned int row = 0; row < mi1.natoms(); row++) {
    unsigned int col=p.my_assig[row];
    if(col<mi2.natoms() && row < mi1.natoms()) {
      v=m[row][col];
      result.mapping1to2[row] = col;
      result.atom_scores[row] = v/10000.0d;
      final_score+=v;
    }
  }
  //printf("Matching Score: %7.3f\n", final_score);
  hungarian_free(&p);
  free_array_to_matrix(m, mi1.natoms());
  free(atom_matrix);
  free(aba_match_scores);

  result.similarity = final_score/10000.0d;
}

NAMSResult * getNAMSResult(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2, const NAMSParameters & params) {
  NAMSResult * result = new NAMSResult;

  result->self_similarity1 = calcSelfSimilarity(molinfo1, params);
  result->self_similarity2 = calcSelfSimilarity(molinfo2, params);
  nams_runner(molinfo1, molinfo2, params, *result, nullptr);
  result->jaccard = result->similarity / ( result->self_similarity1 + result->self_similarity2 - result->similarity );

  return result;
}

double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2, const NAMSParameters & params) {
  std::unique_ptr< NAMSResult > result( getNAMSResult(molinfo1, molinfo2, params) );
  return result->jaccard;
}


} // namespace NAMS
} // namespace RDKit
