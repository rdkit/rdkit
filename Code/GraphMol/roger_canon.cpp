//
//  Copyright (C) 2014 Greg Landrum
//  Adapted from pseudo-code from Roger Sayle
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "roger_canon.h"

#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include <iostream>
#include <cassert>
#include <GraphMol/RankAtoms.h>

namespace RDKit {
  namespace Canon{
    void CreateSinglePartition(unsigned int nAtoms,
                               int *order,
                               int *count,
                               canon_atom *atoms) {
      for( unsigned int i=0; i<nAtoms; i++ ) {
        atoms[i].index = 0;
        order[i] = i;
        count[i] = 0;
      }
      count[0] = nAtoms;
    }

    void ActivatePartitions(unsigned int nAtoms,
                            int *order,
                            int *count,
                            int &activeset,int *next) {
      unsigned int i,j;
      activeset = -1;
      for( i=0; i<nAtoms; i++ )
        next[i] = -2;

      i = 0;
      do {
        j = order[i];
        if( count[j] > 1 ) {
          next[j] = activeset;
          activeset = j;
          i += count[j];
        } else i++;
      } while( i < nAtoms );
    }

    void RankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                      bool includeChirality,bool includeIsotopes) {
      std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        atoms[i].atom=mol.getAtomWithIdx(i);
        atoms[i].index=i;
      }
      AtomCompareFunctor ftor(&atoms.front(),mol);

      canon_atom *data=&atoms.front();
      int *count=(int *)malloc(atoms.size()*sizeof(int));
      int *order=(int *)malloc(atoms.size()*sizeof(int));
      int activeset;
      int *next=(int *)malloc(atoms.size()*sizeof(int));
      CreateSinglePartition(atoms.size(),order,count,data);
      ActivatePartitions(atoms.size(),order,count,activeset,next);
      RefinePartitions(mol,data,ftor,false,order,count,activeset,next);

      ftor.df_useNbrs=true;
      ActivatePartitions(atoms.size(),order,count,activeset,next);
      RefinePartitions(mol,data,ftor,true,order,count,activeset,next);
      BreakTies(mol,data,ftor,true,order,count,activeset,next);

      res.resize(mol.getNumAtoms());
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        res[order[i]]=i;
      }
      free(count); free(order); free(next);

    } // end of RankAtoms()
    
  } // end of Canon namespace
} // end of RDKit namespace
