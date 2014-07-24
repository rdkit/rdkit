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

namespace RDKit {
  namespace Canon{
    void CreateSinglePartition(unsigned int nAtoms,
                               unsigned int *order,
                               unsigned int *count,
                               unsigned int *atomIndices) {
      for( unsigned int i=0; i<nAtoms; i++ ) {
        atomIndices[i] = 0;
        order[i] = i;
        count[i] = 0;
      }
      count[0] = nAtoms;
    }

    void ActivatePartitions(unsigned int nAtoms,
                            unsigned int *order,
                            unsigned int *count,
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

#if 0
    void RefinePartitions(unsigned int nAtoms,
                          func compar, int mode,
                          unsigned int *order,
                          unsigned int *count,
                          int &activeset, int *next,
                          unsigned int *atomIndices)
    {
      register int partition;
      register int symclass;
      register int *start;
      register int offset;
      register int index;
      register int len;
      register int i;

      while( activeset != -1 ) {
        partition = activeset;
        activeset = next[partition];
        next[partition] = -2;

        len = count[partition]; 
        offset = atom[partition].index;
        start = order+offset;
        hanoisort(start,len,count,compar);

        index = start[0];
        for( i=count[index]; i<len; i++ ) {
          index = start[i];
          if( count[index] )
            symclass = offset+i;
          atom[index].index = symclass;

          if( mode ) {
            foreach nbor of atom[index] {
              offset = atom[nbor].index;
              partition = order[offset];
              if( (count[partition]>1) &&
                  (next[partition]==-2) ) {
                next[partition] = activeset;
                activeset = partition;
              }
            }
          }
        }
      }
    }

    void SymmetryPerception( Molecule mol )
    {
      CreateSinglePartition();
      ActivatePartitions();
      RefinePartitions(compare1,false);
      ActivatePartitions();
      RefinePartitions(compare2,true);
    }
#endif

    
  } // end of Canon namespace
} // end of RDKit namespace
