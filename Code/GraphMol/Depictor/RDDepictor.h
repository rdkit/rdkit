//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
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

#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <boost/smart_ptr.hpp>

namespace RDKit {
class ROMol;
}

namespace RDDepict {

RDKIT_DEPICTOR_EXPORT extern bool
    preferCoordGen;  // Ignored if coordgen support isn't active

typedef boost::shared_array<double> DOUBLE_SMART_PTR;

class RDKIT_DEPICTOR_EXPORT DepictException : public std::exception {
 public:
  DepictException(const char *msg) : _msg(msg) {}
  DepictException(const std::string msg) : _msg(msg) {}
  const char *what() const noexcept override { return _msg.c_str(); }
  ~DepictException() noexcept override = default;

 private:
  std::string _msg;
};

//! \brief Set the path to the file containing the ring system templates
/*!

  \param templatePath the file path to a file containing the ring system
  templates. Each template must be a single line in the file represented using
  CXSMILES, and the structure should be a single ring system.

  \throws DepictException if any of the templates are invalid
*/
void RDKIT_DEPICTOR_EXPORT
setRingSystemTemplates(const std::string templatePath);

//! \brief Add ring system templates to be used in 2D coordinater generation.
/// If there are duplicates, the most recently added template will be used.
/*!

  \param templatePath the file path to a file containing the ring system
  templates. Each template must be a single line in the file represented using
  CXSMILES, and the structure should be a single ring system.

  \throws DepictException if any of the templates are invalid
*/
void RDKIT_DEPICTOR_EXPORT
addRingSystemTemplates(const std::string templatePath);

//! \brief Load default ring system templates to be used in 2D coordinate
//! generation
void RDKIT_DEPICTOR_EXPORT loadDefaultRingSystemTemplates();

struct RDKIT_DEPICTOR_EXPORT Compute2DCoordParameters {
  const RDGeom::INT_POINT2D_MAP *coordMap =
      nullptr;  //!< a map of int to Point2D, between atom IDs and their
                //!< locations.  This is the container the user needs to
                //!< fill if he/she wants to specify coordinates for a portion
                //!< of the molecule, defaults to 0
  bool canonOrient = false;  //!< canonicalize the orientation so that the long
                             //!< axes align with the x-axis etc.
  bool clearConfs = true;  //!< clear all existing conformations on the molecule
                           //!< before adding the 2D coordinates instead of
                           //!< simply adding to the list
  unsigned int nFlipsPerSample = 0;  //!< the number of rotatable bonds that are
                                     //!< flipped at random for each sample
  unsigned int nSamples = 0;         //!< the number of samples
  int sampleSeed = 0;                //!< seed for the random sampling process
  bool permuteDeg4Nodes = false;  //!< try permuting the drawing order of bonds
                                  //!< around atoms with four neighbors in order
                                  //!< to improve the depiction
  bool forceRDKit = false;        //!< use RDKit to generate coordinates even if
                                  //!< preferCoordGen is set to true
  bool useRingTemplates = false;  //!< whether to use ring system templates for
                                  //!< generating initial coordinates

  Compute2DCoordParameters() = default;
};

//! \brief Generate 2D coordinates (a depiction) for a molecule
/*!

  \param mol the molecule were are interested in

  \param params parameters used for 2D coordinate generation

  \return ID of the conformation added to the molecule containing the
  2D coordinates

*/
RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoords(
    RDKit::ROMol &mol, const Compute2DCoordParameters &params);

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

  \param forceRDKit - use RDKit to generate coordinates even if
        preferCoordGen is set to true

  \param useRingTemplates whether to use ring system templates for generating
      initial coordinates

  \return ID of the conformation added to the molecule containing the
  2D coordinates

*/
RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoords(
    RDKit::ROMol &mol, const RDGeom::INT_POINT2D_MAP *coordMap = nullptr,
    bool canonOrient = false, bool clearConfs = true,
    unsigned int nFlipsPerSample = 0, unsigned int nSamples = 0,
    int sampleSeed = 0, bool permuteDeg4Nodes = false, bool forceRDKit = false,
    bool useRingTemplates = false);

//! \brief Compute the 2D coordinates such the interatom distances
///  mimic those in a distance matrix
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

  \param forceRDKit - use RDKit to generate coordinates even if
        preferCoordGen is set to true

  \return ID of the conformation added to the molecule containing the
  2D coordinates


*/
RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoordsMimicDistMat(
    RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat = nullptr,
    bool canonOrient = true, bool clearConfs = true, double weightDistMat = 0.5,
    unsigned int nFlipsPerSample = 3, unsigned int nSamples = 100,
    int sampleSeed = 25, bool permuteDeg4Nodes = true, bool forceRDKit = false);

//! \brief Compute 2D coordinates where a piece of the molecule is
///  constrained to have the same coordinates as a reference.
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
  \param acceptFailure - (optional) if true, standard depictions will be
                         generated for molecules that don't have a substructure
                         match to the reference; if false, throws a
                         DepictException.
  \param forceRDKit - (optional) use RDKit to generate coordinates even if
                      preferCoordGen is set to true
  \param allowOptionalAttachments -  (optional) if true, terminal dummy atoms in
                         the reference are ignored if they match an implicit
                         hydrogen in the molecule, and a constrained
                         depiction is still attempted
  RETURNS:

  \return MatchVectType with (queryAtomidx, molAtomIdx) pairs used for
          the constrained depiction
*/
RDKIT_DEPICTOR_EXPORT RDKit::MatchVectType generateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId = -1,
    const RDKit::ROMol *referencePattern =
        static_cast<const RDKit::ROMol *>(nullptr),
    bool acceptFailure = false, bool forceRDKit = false,
    bool allowOptionalAttachments = false);

//! \brief Compute 2D coordinates where a piece of the molecule is
///  constrained to have the same coordinates as a reference.
/*!
  This function generates a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  This is useful for, for example, generating depictions of SAR data
  sets so that the cores of the molecules are all oriented the same way.
  This overload allow to specify the (referenceAtom, molAtom) index pairs
  which should be matched as MatchVectType. Please note that the
  vector can be shorter than the number of atoms in the reference.

  ARGUMENTS:

  \param mol -    the molecule to be aligned, this will come back
                  with a single conformer.
  \param reference -    a molecule with the reference atoms to align to;
                        this should have a depiction.
  \param refMatchVect -  a MatchVectType that will be used to
                         generate the atom mapping between the molecule
                         and the reference.
  \param confId -       (optional) the id of the reference conformation to use
  \param forceRDKit - (optional) use RDKit to generate coordinates even if
                      preferCoordGen is set to true
*/
RDKIT_DEPICTOR_EXPORT void generateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const RDKit::MatchVectType &refMatchVect, int confId = -1,
    bool forceRDKit = false);

//! \brief Generate a 2D depiction for a molecule where all or part of
///  it mimics the coordinates of a 3D reference structure.
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
  \param forceRDKit - (optional) use RDKit to generate coordinates even if
                      preferCoordGen is set to true
*/
RDKIT_DEPICTOR_EXPORT void generateDepictionMatching3DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId = -1,
    RDKit::ROMol *referencePattern = nullptr, bool acceptFailure = false,
    bool forceRDKit = false);

//! \brief Rotate the 2D depiction such that the majority of bonds have an
//! angle with the X axis which is a multiple of 30 degrees.
/*!

  ARGUMENTS:
  \param mol - the molecule to be rotated
  \param confId - (optional) the id of the reference conformation to use
  \param minimizeRotation - (optional) if false (the default), the molecule
  is rotated such that the majority of bonds have an angle with the
  X axis of 30 or 90 degrees. If true, the minimum rotation is applied
  such that the majority of bonds have an angle with the X axis of
  0, 30, 60, or 90 degrees, with the goal of altering the initial
  orientation as little as possible .
*/

RDKIT_DEPICTOR_EXPORT void straightenDepiction(RDKit::ROMol &mol,
                                               int confId = -1,
                                               bool minimizeRotation = false);

//! \brief Normalizes the 2D depiction.
/*!
  If canonicalize is != 0, the depiction is subjected to a canonical
  transformation such that its main axis is aligned along the X axis
  (canonicalize >0, the default) or the Y axis (canonicalize <0).
  If canonicalize is 0, no canonicalization takes place.
  If scaleFactor is <0.0 (the default) the depiction is scaled such
  that bond lengths conform to RDKit standards. The applied scaling
  factor is returned.

  ARGUMENTS:
  \param mol          - the molecule to be normalized
  \param confId       - (optional) the id of the reference conformation to use
  \param canonicalize - (optional) if != 0, a canonical transformation is
                        applied: if >0 (the default), the main molecule axis is
                        aligned to the X axis, if <0 to the Y axis.
                        If 0, no canonical transformation is applied.
  \param scaleFactor  - (optional) if >0.0, the scaling factor to apply. The
                        default (-1.0) means that the depiction is automatically
                        scaled such that bond lengths are the standard RDKit
                        ones.
  RETURNS:

  \return the applied scaling factor.
*/

RDKIT_DEPICTOR_EXPORT double normalizeDepiction(RDKit::ROMol &mol,
                                                int confId = -1,
                                                int canonicalize = 1,
                                                double scaleFactor = -1.0);
};  // namespace RDDepict

#endif
