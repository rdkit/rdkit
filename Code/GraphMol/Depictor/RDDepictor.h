//
//  Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDDEPICTOR_H
#define RDDEPICTOR_H

#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <boost/smart_ptr.hpp>

namespace RDKit {
class ROMol;
}

namespace RDDepict {

#ifdef RDK_BUILD_COORDGEN_SUPPORT
RDKIT_DEPICTOR_EXPORT extern bool preferCoordGen;
#endif

typedef boost::shared_array<double> DOUBLE_SMART_PTR;

class RDKIT_DEPICTOR_EXPORT DepictException : public std::exception {
 public:
  DepictException(const char *msg) : _msg(msg){};
  DepictException(const std::string msg) : _msg(msg){};
  const char *what() const noexcept override { return _msg.c_str(); };
  ~DepictException() noexcept {};

 private:
  std::string _msg;
};

//! \brief Generate 2D coordinates (a depiction) for a molecule
/*!

  \param mol the molecule were are interested in

  \param coordMap a map of int to Point2D, between atom IDs and
  their locations.  This is the container the user needs to fill if
  he/she wants to specify coordinates for a portion of the molecule,
  defaults to 0

  \param canonOrient canonicalize the orientation so that the long
  axes align with the x-axis etc.

  \param clearConfs clear all existing conformations on the molecule
  before adding the 2D coordinates instead of simply adding to the
  list

  \param nFlipsPerSample - the number of rotatable bonds that are
  flipped at random for each sample

  \param nSamples - the number of samples

  \param sampleSeed - seed for the random sampling process

  \param permuteDeg4Nodes - try permuting the drawing order of bonds around
        atoms with four neighbors in order to improve the depiction

  \return ID of the conformation added to the molecule containing the
  2D coordinates

*/
RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoords(
    RDKit::ROMol &mol, const RDGeom::INT_POINT2D_MAP *coordMap = nullptr,
    bool canonOrient = false, bool clearConfs = true,
    unsigned int nFlipsPerSample = 0, unsigned int nSamples = 0,
    int sampleSeed = 0, bool permuteDeg4Nodes = false, bool forceRDKit = false);

//! \brief Compute the 2D coordinates such the interatom distances
//   mimic those in a distance matrix
/*!

  This function generates 2D coordinates such that the inter-atom
  distances mimic those specified via dmat. This is done by randomly
  sampling(flipping) the rotatable bonds in the molecule and
  evaluating a cost function which contains two components. The
  first component is the sum of inverse of the squared inter-atom
  distances, this helps in spreading the atoms far from each
  other. The second component is the sum of squares of the
  difference in distance between those in dmat and the generated
  structure.  The user can adjust the relative importance of the two
  components via a adjustable parameter (see below)

  ARGUMENTS:

  \param mol - molecule to generate coordinates for

  \param dmat - the distance matrix we want to mimic, this is a
  symmetric N by N matrix where N is the number of atoms in mol. All
  negative entries in dmat are ignored.

  \param canonOrient - canonicalize the orientation after the 2D
  embedding is done

  \param clearConfs - clear any previously existing conformations on
  mol before adding a conformation

  \param weightDistMat - A value between 0.0 and 1.0, this
  determines the importance of mimicing the inter atoms
  distances in dmat. (1.0 - weightDistMat) is the weight associated
  to spreading out the structure (density) in the cost function

  \param nFlipsPerSample - the number of rotatable bonds that are
  flipped at random for each sample

  \param nSamples - the number of samples

  \param sampleSeed - seed for the random sampling process

  \param permuteDeg4Nodes - try permuting the drawing order of bonds around
        atoms with four neighbors in order to improve the depiction

  \return ID of the conformation added to the molecule containing the
  2D coordinates


*/
RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoordsMimicDistMat(
    RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat = nullptr,
    bool canonOrient = true, bool clearConfs = true, double weightDistMat = 0.5,
    unsigned int nFlipsPerSample = 3, unsigned int nSamples = 100,
    int sampleSeed = 25, bool permuteDeg4Nodes = true, bool forceRDKit = false);

//! \brief Compute 2D coordinates where a piece of the molecule is
//   constrained to have the same coordinates as a reference.
/*!
  This function generates a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  This is useful for, for example, generating depictions of SAR data
  sets so that the cores of the molecules are all oriented the same way.

  ARGUMENTS:

  \param mol -    the molecule to be aligned, this will come back
                  with a single conformer.
  \param reference -    a molecule with the reference atoms to align to;
                        this should have a depiction.
  \param confId -       (optional) the id of the reference conformation to use
  \param referencePattern -  (optional) a query molecule to be used to
                             generate the atom mapping between the molecule
                             and the reference.
  \param acceptFailure - (optional) if True, standard depictions will be
  generated
                         for molecules that don't have a substructure match to
  the
                         reference; if false, throws a DepictException.

*/
RDKIT_DEPICTOR_EXPORT void generateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId = -1,
    RDKit::ROMol *referencePattern = static_cast<RDKit::ROMol *>(nullptr),
    bool acceptFailure = false, bool forceRDKit = false);

//! \brief Generate a 2D depiction for a molecule where all or part of
//   it mimics the coordinates of a 3D reference structure.
/*!
  Generates a depiction for a molecule where a piece of the molecule
  is constrained to have coordinates similar to those of a 3D reference
  structure.

  ARGUMENTS:
  \param mol - the molecule to be aligned, this will come back
               with a single conformer containing 2D coordinates
  \param reference - a molecule with the reference atoms to align to.
                     By default this should be the same as mol, but with
                     3D coordinates
  \param confId - (optional) the id of the reference conformation to use
  \param refPattern - (optional) a query molecule to map a subset of
                      the reference onto the mol, so that only some of the
                      atoms are aligned.
  \param acceptFailure - (optional) if true, standard depictions will be
  generated
                         for molecules that don't match the reference or the
                         referencePattern; if false, throws a DepictException.
*/
RDKIT_DEPICTOR_EXPORT void generateDepictionMatching3DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId = -1,
    RDKit::ROMol *referencePattern = nullptr, bool acceptFailure = false,
    bool forceRDKit = false);
};  // namespace RDDepict

#endif
