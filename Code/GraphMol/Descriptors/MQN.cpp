// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <boost/foreach.hpp>
#include <vector>
#include <algorithm>

namespace RDKit{
  namespace Descriptors {
    std::vector<unsigned int> calcMQNs(const ROMol &mol,
                                      bool force){
      // FIX: use force value to enable caching
      std::vector<unsigned int> res(42,0);

      // ---------------------------------------------------
      // atom-centered things
      // Note: We're not doing exactly the same thing
      //       as the original paper on polarity counts
      //       since we're using different donor and acceptor
      //       definitions.
      ROMol::VERTEX_ITER atBegin,atEnd;
      boost::tie(atBegin,atEnd) = mol.getVertices();  
      while(atBegin!=atEnd){
        const ATOM_SPTR at=mol[*atBegin];
        ++atBegin;
        unsigned int nHs = at->getTotalNumHs();
        unsigned int nRings=mol.getRingInfo()->numAtomRings(at->getIdx());
        switch(at->getAtomicNum()){
        case 0:
        case 1:
          break;
        case 6:
          res[0]++;break;
        case 9:
          res[1]++;break;
        case 17:
          res[2]++;break;
        case 35:
          res[3]++;break;
        case 53:
          res[4]++;break;
        case 16:
          res[5]++;break;
        case 15:
          res[6]++;break;
        case 7:
          if(!nRings) res[7]++;
          else res[8]++;
          if(at->getDegree()!=4){
            res[19]++; // number of acceptor sites
            res[20]++; // number of acceptor atoms
          }
          if(nHs){
            res[21]+=nHs; // number of donor sites
            res[22]++; // number of donor atoms
          }
          break;
        case 8:
          if(!nRings) res[9]++;
          else res[10]++;
          res[20]++; // number of acceptor atoms
          if(at->getFormalCharge()!=-1){
            res[19]+=2; // number of acceptor sites
          } else {
            res[19]+=3; // number of acceptor sites
          }
          if(nHs){
            res[21]+=nHs; // number of donor sites
            res[22]++; // number of donor atoms
          }
          break;
        default:
          res[11]++;break;
        }

        if(at->getFormalCharge()>0) res[24]++; // positive charges
        else if(at->getFormalCharge()<0) res[23]++; // negative charges

        if(at->getAtomicNum()!=1){
          switch(at->getDegree()){
          case 1:
            res[25]++;
            break;
          case 2:
            if(!nRings) res[26]++;
            else res[29]++;
            break;
          case 3:
            if(!nRings) res[27]++;
            else res[30]++;
            break;
          case 4:
            if(!nRings) res[28]++;
            else res[31]++;
            break;
          }
          if(nRings>=2) res[40]++;
        }
      }
      
      // ---------------------------------------------------
      // bond counts:
      unsigned int nAromatic=0;
      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
        const BOND_SPTR bond = mol[*firstB];
        if(bond->getIsAromatic()) ++nAromatic;
        unsigned int nRings=mol.getRingInfo()->numBondRings(bond->getIdx());
        switch(bond->getBondType()){
        case Bond::SINGLE:
          if(!nRings) res[12]++;
          else res[15]++;
          break;
        case Bond::DOUBLE:
          if(!nRings) res[13]++;
          else res[16]++;
          break;
        case Bond::TRIPLE:
          if(!nRings) res[14]++;
          else res[17]++;
          break;
        default:
          break;
        }
        if(nRings>=2) res[41]++;
        ++firstB;
      }
      // rather than do the work to kekulize the molecule, we cheat
      // by just dividing the number of aromatic bonds evenly among the
      // cyclic single bond and cyclic double bond bins and give any
      // remainder to the single bonds
      res[15] += nAromatic/2;
      res[16] += nAromatic/2;
      if(nAromatic%2) res[15]++;
      res[18] = calcNumRotatableBonds(mol);

      // ---------------------------------------------------
      //  ring size counts
      BOOST_FOREACH(const INT_VECT &iv,mol.getRingInfo()->atomRings()){
        if(iv.size()<10){
          res[iv.size()+29]++;
        } else {
          res[39]++;
        }
      }
      
      return res;
    }
  } // end of namespace Descriptors
} // end of namespace RDKit
