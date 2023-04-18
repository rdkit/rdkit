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
// (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

#include <GraphMol/Fingerprints/NAMS.h>
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
  nring1 = mol.getRingInfo()->numAtomRings(at1);
  nring2 = mol.getRingInfo()->numAtomRings(at2);
  if ( do_isomerism ) { /*TODO: implement chirality?*/ }
  const Bond* bond = mol.getBondBetweenAtoms(at1,at2);
  assert( bond != nullptr );
  ringbond = mol.getRingInfo()->minBondRingSize(bond->getIdx()) != 0;
  aromatic = bond->getIsAromatic();
  order = bond->getBondTypeAsDouble();
  if ( do_isomerism ) { /*TODO: implement isomerism?*/ }
}

bool
BondInfoType::operator<( const BondInfoType & r ) const {
  return std::make_tuple(ele1, nring1, chr1, ele2, nring2, chr2, ringbond, aromatic, order, dbstero12)
    < std::make_tuple(r.ele1, r.nring1, r.chr1, r.ele2, r.nring2, r.chr2, r.ringbond, r.aromatic, r.order, r.dbstero12);

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
 Returns a pre-processed data structure for the molecule
 Top level is a vector with one entry per atom
 Second level is a vector with one entry per distance step
 Third level is a vector of bond information for all the bonds at that level
*/
std::vector< std::vector< std::vector< BondInfo > > >
get_connection_matrix(const ROMol &mol, bool do_isomerism = false) {
  // From nams-docker/makenamsdb/recoder.py:Recoder::get_mol_info()

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

std::ostream & operator<<( std::ostream & os, BondInfoType const & type ) {
  os << type.ele1 << ' ' << type.nring1 << ' ' << type.chr1;
  os << ' ' << type.ele2 << ' ' << type.nring2 << ' ' << type.chr2;
  os << ' ' << type.ringbond << ' ' << type.aromatic << ' ' << type.order*10 << ' ' << type.dbstero12;
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
  ss << mat_aba_types.size() << ' ' << nbonds << ' ' << aba_types.size() << '\n';
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

double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2) {
  // For debugging right now
  std::cerr << "--%%----------------------------------\n";
  std::cerr << molinfo1.dump(1);
  std::cerr << molinfo2.dump(2);
  std::cerr << "--%%----------------------------------\n";
  return 0;
}

} // namespace NAMS
} // namespace RDKit
