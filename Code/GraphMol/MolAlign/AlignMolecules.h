//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_ALIGNMOLECULES_H_
#define _RD_ALIGNMOLECULES_H_

#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <vector>

namespace RDKit {
  typedef std::vector< std::pair<int,int> > MatchVectType; 

  class Conformer;
  class ROMol;
  namespace MolAlign {
    class MolAlignException : public std::exception {
    public:
      //! construct with an error message
      MolAlignException(const char *msg) : _msg(msg) {};
      //! construct with an error message
      MolAlignException(const std::string msg) : _msg(msg) {};
      //! get the error message
      const char *message () const { return _msg.c_str(); };
      ~MolAlignException () throw () {};
    private:
      std::string _msg;
    };

    //! Alignment functions

    //! Compute the transformation required to align a molecule
    /*!
      The 3D transformation required to align the specied conformation in the probe molecule
      to a specified conformation in the reference molecule is computed so that the root mean
      squared distance between a specified set of atoms is minimized

      \param prbMol    molecule that is to be aligned
      \param refMol    molecule used as the reference for the alignment
      \param trans     storage for the computed transform
      \param prbCid    ID of the conformation in the probe to be used 
                       for the alignment (defaults to first conformation)
      \param refCid    ID of the conformation in the ref molecule to which 
                       the alignment is computed (defaults to first conformation)
      \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                       used to compute the alignments. If this mapping is 
                       not specified an attempt is made to generate on by
                       substructure matching
      \param weights   Optionally specify weights for each of the atom pairs
      \param reflect   if true reflect the conformation of the probe molecule
      \param maxIters  maximum number of iteration used in mimizing the RMSD

      <b>Returns</b>
      RMSD value
    */
    double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol, 
                                         RDGeom::Transform3D &trans,
                                         int prbCid=-1, int refCid=-1, 
                                         const MatchVectType *atomMap=0, 
                                         const RDNumeric::DoubleVector *weights=0, 
                                         bool reflect=false,
                                         unsigned int maxIters=50);
    
    
    //! Optimally (minimum RMSD) align a molecule to another molecule
    /*!
      The 3D transformation required to align the specied conformation in the probe molecule
      to a specified conformation in the reference molecule is computed so that the root mean
      squared distance between a specified set of atoms is minimized. This transforms is them
      applied to the specified conformation in the probe molecule
      
      \param prbMol    molecule that is to be aligned
      \param refMol    molecule used as the reference for the alignment
      \param prbCid    ID of the conformation in the probe to be used 
                       for the alignment (defaults to first conformation)
      \param refCid    ID of the conformation in the ref molecule to which 
                       the alignment is computed (defaults to first conformation)
      \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                       used to compute the alignments. If this mapping is 
                       not specified an attempt is made to generate on by
                       substructure matching
      \param weights   Optionally specify weights for each of the atom pairs
      \param reflect   if true reflect the conformation of the probe molecule
      \param maxIters  maximum number of iteration used in mimizing the RMSD

      <b>Returns</b>
      RMSD value
    */
    double alignMol(ROMol &prbMol, const ROMol &refMol, 
                            int prbCid=-1, int refCid=-1,
                            const MatchVectType *atomMap=0, 
                            const RDNumeric::DoubleVector *weights=0, 
                            bool reflect=false, unsigned int maxIters=50);

    //! \brief Calculates eigenvalues (= principal moments of inertia) and eigenvectors of inertia tensor
    /*!
      \param mol                molecule that is to be used
      \param eigenVals          array[3] for storing eigenvalues of inertia tensor
      \param eigenVects         array[3][3] for storing eigenVecs of inertia tensor
      \param confId             ID of the conformation in the probe to be used
                                for the alignment (defaults to first conformation)
      \param weights            Ptr to vector of weights of atom positions to be used.
                                                defaults to masses of atoms (if weights==NULL Ptr)
      \param maxIterations      maximum number of iterations used to align molecule
     */
    void getMomentsOfInertia(const ROMol &mol, double eigenVals[3], double eigenVecs[3][3], int confId = -1,
                             RDNumeric::DoubleVector *weights=0, unsigned int maxIterations=50);

    //! \brief Calculate alignment transform necessary for aligning a molecule to its principal moments of inertia
    /*!
      \param mol		molecule that is to be aligned
      \param trans	    	a RDGeom::Transform3D object (affine transformation matrix)
      \param eigenVals          an array of doubles to store the eigenvalues in (= principal moments of inertia)
                                defaults to NULL (nothing is stored)
      \param eigenVecs          a 3x3 array of doubles to store the eigenvectors in
                                defaults to NULL (nothing is stored)
      \param confId		ID of the conformation in the probe to be used
                       		for the alignment (defaults to first conformation)
      \param weights		Ptr to vector of weights of atom positions to be used.
      	  	  	  	  		defaults to masses of atoms (if weights==NULL Ptr)
      \param maxIterations 	maximum number of iterations used to align molecule
     */
    void getPrincAxesTransform(ROMol &mol, RDGeom::Transform3D &trans, double eigenVals[3]=NULL, double eigenVecs[3][3]=NULL,
                               int confId = -1, RDNumeric::DoubleVector *weights=0, unsigned int maxIterations=50);

    //! Align the conformations of a molecule using a common set of atoms
    /*! 
      \param mol       The molecule of interest
      \param atomIds   vector of atoms to be used to generate the alignment.
                       All atoms will be used is not specified
      \param confIds   vector of conformations to align - defaults to all
      \param weights   vector of weights to applied to particular atom pairs
                       defaults to all weights = 1
      \param reflect   toggles reflecting (about the origin) the alignment
      \param maxIters  the maximum number of iterations to attempt
      \param RMSlist   if nonzero, this will be used to return the RMS values
                       between the reference conformation and the other aligned
                       conformations
    */
    void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds=0,
                            const std::vector<unsigned int> *confIds=0,
                            const RDNumeric::DoubleVector *weights=0, 
                            bool reflect=false, unsigned int maxIters=50,
                            std::vector<double> *RMSlist=0);
  }
}
#endif
