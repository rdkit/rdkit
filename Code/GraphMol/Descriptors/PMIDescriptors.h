/*
 * PMIDescriptors.h
 *
 *  Created on: Jun 25, 2014
 *      Author: hahnda6
 */

#ifndef PMIDESCRIPTORS_H_
#define PMIDESCRIPTORS_H_

#include <GraphMol/ROMol.h>
#include <Numerics/Vector.h>

namespace RDKit{
  namespace Descriptors{

    struct Moments {
      double PMI1, PMI2, PMI3, NPR1, NPR2;
      inline double &operator[] (unsigned int i) {
        switch (i) {
          case 0:
            return PMI1;
            break;
          case 1:
            return PMI2;
            break;
          case 2:
            return PMI3;
            break;
          case 3:
            return NPR1;
            break;
          case 4:
            return NPR2;
            break;
          default:
            throw IndexErrorException(i);
        }
      }
    };
    //! \brief Calculate PMI descriptors for a molecule
    /*!
    \param mol   	    molecule that is to be used
    \param PMI              a 5x1 array of doubles to store results in (Ix, Iy, Iz, Ix/Iz, Iy/Iz)
    \param confId	    ID of the conformation in the probe to be used
                     	    for the alignment (defaults to first conformation)
    \param weights	    Ptr to vector of weights of atom positions to be used.
    	  	  	    defaults to masses of atoms (if weights==NULL Ptr)
    \param maxIterations    maximum number of iterations used to align molecule
     */
    void calcPMIDescriptors(const RDKit::ROMol &mol, Moments &MoI, int confId=-1, RDNumeric::DoubleVector *weights=0,
                               unsigned int maxIterations=50);

  }
}

#endif /* PMIDESCRIPTORS_H_ */
