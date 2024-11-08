//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDMIF_DESCRIPTORS_H
#define RDMIF_DESCRIPTORS_H

#include <vector>
#include <Geometry/point.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <ForceField/MMFF/Nonbonded.h>
#include <ForceField/UFF/Nonbonded.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/RDKitBase.h>

namespace RDMIF {
//! \brief constructs a UniformRealValueGrid3D which fits to the molecule mol
/*!
 \returns a pointer to a UniformRealValueGrid which has a spacing of \c spacing
 (default = 0.5 Angstrom) and a spatial extent of conformer \c confId
 (default=-1) of molecule \c mol plus a defined margin \c margin on each side
 (default = 5.0 Angstrom).

 \param mol		Molecule for which the grid is constructed
 \param confId		Id of Conformer which is used
 \param margin 	distance from molecule to the grid's wall
 \param spacing	grid's spacing

 <b>Returns</b>
 pointer to a grid
 */
std::unique_ptr<RDGeom::UniformRealValueGrid3D> constructGrid(
    const RDKit::ROMol &mol, int confId = -1, double margin = 5.0,
    double spacing = 0.5);

//! \brief calculates a descriptor at every grid point of MIF
/*!
 any unary function/functor can be used taking a grid point of type
 RDKit::Point3D as parameter \param grd			a UniformRealValueGrid3D
 \param functor		any unary function/functor which takes four doubles as
 parameters (3 coordinates, 1 threshold) \param thres	                cuts of
 atom contributions if distance of interaction higher than threshold, negative
 threshold: no threshold, defaults to -1.0
 */
template <typename T>
void calculateDescriptors(RDGeom::UniformRealValueGrid3D &grd, T functor,
                          double thres = -1.0) {
  const RDGeom::Point3D &offSet = grd.getOffset();
  double oX = offSet.x, oY = offSet.y, z = offSet.z;
  double x = oX, y = oY;
  double spacing = grd.getSpacing();
  if (thres < 0) {
    thres = spacing * grd.getSize();
  }
  thres *= thres;  // comparing squared distance

  unsigned int id = 0;
  const boost::shared_array<double> &data = grd.getDataPtr();

  for (unsigned int idZ = 0; idZ < grd.getNumZ(); idZ++) {
    for (unsigned int idY = 0; idY < grd.getNumY(); idY++) {
      for (unsigned int idX = 0; idX < grd.getNumX(); idX++) {
        data[id++] = functor(x, y, z, thres);
        x += spacing;
      }
      y += spacing;
      x = oX;
    }
    z += spacing;
    y = oY;
  }
}

//! \brief class for calculation of distance to closest atom of \c mol from grid
//! point \c pt
class DistanceToClosestAtom {
 public:
  //! construct a distanceToClosestAtom from a ROMol
  DistanceToClosestAtom(const RDKit::ROMol &mol, int confId = -1);

  //! Calculates the closest distance, \returns distance in Angstrom
  double operator()(double x, double y, double z, double thres);

 private:
  unsigned int d_nAtoms;
  std::vector<double> d_pos;
};

//! \brief class for calculation of electrostatic interaction (Coulomb energy)
//! between probe and molecule in vaccuum (no dielectric)
/*
 d_nAtoms     No of atoms in molecule
 d_softcore   true if softcore interaction is used, else minimum cutoff distance
 is used d_probe      charge of probe [e] d_absVal     if true, absolute values
 of interactions are calculated d_alpha      softcore interaction parameter
 [A^2] d_cutoff     squared minimum cutoff distance [A^2] d_charges    vector of
 doubles with all partial charges of the atoms in a molecule [e] d_pos vector of
 Point3Ds with all positions of the atoms in a molecule [A] d_prefactor
 prefactor taking into account the geometric factors,natural constants and
 conversion of units
 */
class Coulomb {
 public:
  Coulomb() = default;
  ~Coulomb() = default;

  //! \brief constructs Coulomb object from vectors of charges and positions
  /*!
   \param charges     vector of charges [e]
   \param pos         vector of postions [A]
   \param probecharge charge of probe [e] (default: 1.0)
   \param absVal      if true, negative (favored) values of interactions are
   calculated (default: false)
   \param alpha       softcore interaction parameter
   [A^2], if zero, a minimum cutoff distance is used (default=0.0)
   \param cutoff
   minimum cutoff distance [A] (default:1.0)

   */
  Coulomb(const std::vector<double> &charges,
          const std::vector<RDGeom::Point3D> &positions,
          double probecharge = 1.0, bool absVal = false, double alpha = 0.0,
          double cutoff = 1.0);

  //! \brief constructs Coulomb object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.

   \param mol         molecule object
   \param confId      conformation id which is used to get positions of atoms
   (default=-1)
   \param absVal      if true, absolute values of interactions are
   calculated (default: false)
   \param probecharge charge of probe [e]
   (default: 1.0)
   \param prop	     property key for retrieving partial charges
   of atoms (default="_GasteigerCharge")
   \param alpha       softcore interaction parameter [A^2], if zero, a minimum
   cutoff distance is used (default=0.0)
   \param cutoff      minimum cutoff distance [A] (default:1.0)

   */
  Coulomb(const RDKit::ROMol &mol, int confId = -1, double probecharge = 1.0,
          bool absVal = false, const std::string &prop = "_GasteigerCharge",
          double alpha = 0.0, double cutoff = 1.0);

  //! \brief calculated the electrostatic interaction at point \c pt in the
  //! molecules field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.
   <b>Returns</b>
   electrostatic interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres);

 private:
  unsigned int d_nAtoms = 0;
  bool d_softcore = false;
  bool d_absVal = false;
  double d_cutoff = 0.00001;
  double d_probe = 1.0;
  double d_alpha = 0.0;
  std::vector<double> d_charges;
  std::vector<double> d_pos;
};

//! \brief class for calculation of electrostatic interaction (Coulomb energy)
//! between probe and molecule by taking a distance-dependent dielectric into
//! account:
/*!
 Same energy term as used in GRID
 References:
 J. Med. Chem. 1985, 28, 849.
 J. Comp. Chem. 1983, 4, 187.
 */
/*
 member variables:
 d_nAtoms:     No of atoms in molecule
 d_softcore:   true if softcore interaction is used, else minimum cutoff
 distance is used d_dielectric: factor of dielectric constants d_probe: charge
 of probe [e] d_absVal:     if true, negative (favored) values of interactions
 are calculated (default: false) d_epsilon:    relative permittivity of solvent
 d_xi:         relative permittivity of solute
 d_alpha:      softcore interaction parameter [A^2]
 d_cutoff:     squared minimum cutoff distance [A^2]
 d_charges:    vector of doubles with all partial charges of the atoms in a
 molecule [e] d_pos:        vector of Point3Ds with all positions of the atoms
 in a molecule [A] d_prefactor:  prefactor taking into account the geometric
 factors,natural constants and conversion of units
 */
class CoulombDielectric {
 public:
  CoulombDielectric()
      : d_nAtoms(0),
        d_softcore(false),
        d_absVal(false),
        d_cutoff(0.001),
        d_probe(1),
        d_epsilon(1.0),
        d_xi(1.0),
        d_alpha(0.0) {
    d_dielectric = (d_xi - d_epsilon) / (d_xi + d_epsilon);
  }

  //! \brief constructs CoulombDielectric object from vectors of charges and
  //! positions
  /*!
   \param charges     vector of charges [e]
   \param pos         vector of postions [A]
   \param probecharge charge of probe [e] (default: 1.0)
   \param absVal      if true, negative (favored) values of interactions are
   calculated (default: false) \param alpha       softcore interaction parameter
   [A^2], if zero, a minimum cutoff distance is used (default=0.0) \param cutoff
   minimum cutoff distance [A] (default:1.0) \param epsilon     relative
   permittivity of solvent (default=80.0) \param xi          relative
   permittivity of solute (default=4.0)

   */
  CoulombDielectric(const std::vector<double> &charges,
                    const std::vector<RDGeom::Point3D> &positions,
                    double probecharge = 1.0, bool absVal = false,
                    double alpha = 0.0, double cutoff = 1.0,
                    double epsilon = 80.0, double xi = 4.0);

  //! \brief constructs Coulomb object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.

   \param mol         molecule object
   \param confId      conformation id which is used to get positions of atoms
   (default=-1) \param probecharge charge of probe [e] (default: 1.0) \param
   absVal      if true, negative (favored) values of interactions are calculated
   (default: false) \param prop	       property key for retrieving partial
   charges of atoms (default="_GasteigerCharge") \param alpha       softcore
   interaction parameter [A^2], if zero, a minimum cutoff distance is used
   (default=0.0) \param cutoff      minimum cutoff distance [A] (default:1.0)
   \param epsilon     relative permittivity of solvent (default=80.0)
   \param xi          relative permittivity of solute (default=4.0)
   */
  CoulombDielectric(const RDKit::ROMol &mol, int confId = -1,
                    double probecharge = 1.0, bool absVal = false,
                    const std::string &prop = "_GasteigerCharge",
                    double alpha = 0.0, double cutoff = 1.0,
                    double epsilon = 80.0, double xi = 4.0);

  //! \brief returns the electrostatic interaction at point \c pt in the
  //! molecule's field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared threshold distance
   <b>Returns</b>
   electrostatic interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres);

 private:
  unsigned int d_nAtoms;
  bool d_softcore, d_absVal;
  double d_cutoff, d_probe;
  double d_epsilon, d_xi, d_alpha;
  double d_dielectric;
  std::vector<double> d_charges, d_sp;
  std::vector<double> d_pos;
  std::vector<double> d_dists;
};

//! \brief class for calculation of Van der Waals interaction between probe and
//! molecule at gridpoint \c pt
/*
 Either the MMFF94 or the UFF VdW term can be calculated with a defined probe
 atom and molecule. d_cutoff:     minimum cutoff distance [A] d_nAtoms:     No
 of atoms in molecule d_R_star_ij   vector of Lennard Jones parameters for every
 interaction (vdW-radii) d_wellDepth   vector of Lennard Jones parameters for
 every interaction (well depths) d_pos	        vector of positions of atoms in
 molecule
 */
class VdWaals {
 public:
  VdWaals() : d_cutoff(1), d_nAtoms(0), d_getEnergy(NULL) {};
  VdWaals(RDKit::ROMol &mol, int confId, unsigned int probeAtomTypeMMFF,
          const std::string &probeAtomTypeUFF, const std::string &FF,
          bool scaling, double cutoff);
  ~VdWaals() {};

  //! \brief returns the VdW interaction at point \c pt in the molecules field
  //! in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.
   <b>Returns</b>
   vdW interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres);

 private:
  double d_cutoff;
  unsigned int d_nAtoms;
  std::vector<double> d_R_star_ij, d_wellDepth;
  std::vector<double> d_pos;

  double (*d_getEnergy)(double, double,
                        double);  // pointer to energy function (MMFF94 or UFF)
  static double d_calcUFFEnergy(double, double, double);  // UFF energy function
  static double d_calcMMFFEnergy(double, double,
                                 double);  // MMFF energy function
};

// Factory functions for VdWaals

//! \brief constructs VdWaals object which uses MMFF94 from a molecule object
/*!
 The molecule \c mol needs to have partial charges set as property of atoms.
 \param mol           molecule object
 \param confId        conformation id which is used to get positions of atoms
 (default=-1) \param probeAtomType MMFF94 atom type for the probe atom
 (default=6, sp3 oxygen) \param scaling       scaling of VdW parameters to take
 hydrogen bonds into account (default=false) \param cutoff        minimum cutoff
 distance [A] (default:1.0)
 */
VdWaals constructVdWaalsMMFF(RDKit::ROMol &mol, int confId = -1,
                             unsigned int probeAtomType = 6,
                             bool scaling = false, double cutoff = 1.0);

//! \brief constructs VdWaals object which uses UFF from a molecule object
/*!
 The molecule \c mol needs to have partial charges set as property of atoms.
 \param mol           molecule object
 \param confId        conformation id which is used to get positions of atoms
 (default=-1) \param probeAtomType UFF atom type for the probe atom (eg. "O_3",
 i.e. a tetrahedral oxygen) \param cutoff        minimum cutoff distance [A]
 (default:1.0)
 */
VdWaals constructVdWaalsUFF(
    RDKit::ROMol &mol, int confId = -1,
    const std::string &probeAtomType = "O_3",
    double cutoff = 1.0);  // constructor for UFF LJ interaction

//! \brief class for calculation of hydrogen bond potential between probe and
//! molecule
/*
 implementation of GRID molecular H-Bond descriptor:
 J.Med.Chem. 1989, 32, 1083.
 J.Med.Chem. 1993, 36, 140.
 J.Med.Chem. 1993, 36, 148.
 member variables:
 d_cutoff:       squared minimum cutoff distance [A^2]
 d_nInteract:    No of interactions between probe and molecule
 d_DAprop:		  either 'A' or 'D', defines whether target atoms (in
 molecule) have to be acceptors or donors for interacting with probe
 d_probetype:    defines type of probe atom, either N or O
 d_targettypes:  vectors of types of target atoms (in molecule), either N or O
 d_pos:          vector of positions of target atoms in molecule
 d_direction:    if target accepting, vector of directions of lone pairs on
 target atoms (if there are two lone pairs, the resulting direction is used) if
 target donating, the X-H bond direction is saved d_plane:        vector of
 plane vectors in which the lone pairs are
 */
class HBond {
 public:
  HBond() : d_cutoff(1.0), d_probetype(O), d_nInteract(0) {};

  //! \brief constructs HBond object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.
   \param mol           molecule object
   \param confId        conformation id which is used to get positions of atoms
   (default=-1) \param probeType     atom type for the probe atom (either, "OH",
   "O", "NH" or "N") \param fixed         for some groups, two different angle
   dependencies are defined in GRID: one which takes some flexibility of groups
   (rotation/swapping of lone pairs and hydrogen) into account and one for
   strictly fixed conformations if true, strictly fixed conformations
   (default=fixed) \param cutoff        minimum cutoff distance [A]
   (default:1.0)
   */
  HBond(const RDKit::ROMol &mol, int confId = -1,
        const std::string &probeType = "OH", bool fixed = true,
        double cutoff = 1.0);
  ~HBond() {};

  //! \brief returns the hydrogen bonding interaction at point \c pt in the
  //! molecules field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.
   <b>Returns</b>
   hydrogen bonding interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres);

  unsigned int getNumInteractions() const { return d_nInteract; }

 private:
  boost::uint8_t d_DAprop;

  enum atomtype { N, O };
  double d_cutoff;
  atomtype d_probetype;
  unsigned int d_nInteract;  // number of HBond interactions

  std::vector<atomtype> d_targettypes;
  std::vector<double> d_pos;         // vector of positions of target atoms
  std::vector<double> d_direction;   // vector of sum of lone pair vectors
  std::vector<double> d_lengths;     // vector of lengths of direction vectors
  std::vector<double> d_plane;       // vector of lone pair plane vectors
  std::vector<double> d_eneContrib;  // energy contributions of all interactions
  std::vector<double> d_vectTargetProbe;  // hydrogen bond direction of probe

  // constructing functions
  unsigned int findSpecials(const RDKit::ROMol &mol, int confId, bool fixed,
                            std::vector<unsigned int> &specials);
  unsigned int findAcceptors(const RDKit::ROMol &mol, int confId,
                             const std::vector<unsigned int> &specials);
  unsigned int findDonors(const RDKit::ROMol &mol, int confId,
                          const std::vector<unsigned int> &specials);
  unsigned int findAcceptors_unfixed(const RDKit::ROMol &mol, int confId,
                                     const std::vector<unsigned int> &specials);
  unsigned int findDonors_unfixed(const RDKit::ROMol &mol, int confId,
                                  const std::vector<unsigned int> &specials);

  void addVectElements(atomtype type, double (*funct)(double, double, double),
                       const RDGeom::Point3D &pos, const RDGeom::Point3D &dir,
                       const RDGeom::Point3D &plane = RDGeom::Point3D(0.0, 0.0,
                                                                      0.0));

  void normalize(double &x, double &y, double &z) const;
  double angle(double x1, double y1, double z1, double x2, double y2,
               double z2) const;

  std::vector<double (*)(double, double, double)> d_function;
};

//! \brief class for calculation of hydrophilic field of molecule
/*
 interaction energy of hydrogen and oxygen of water with the molecule is
 calculated at each point as a hydrogen bond interaction (either OH or O probe).
 the favored interaction is returned. member variables: d_hbondO
 HBond descriptor with "O" probe d_hbondOH				   HBond
 descriptor with "OH" probe
 */
class Hydrophilic {
 public:
  //! \brief constructs Hydrophilic object from a molecule object
  /*!
   params:
   \param mol           		 molecule object
   \param confId        		 conformation id which is used to get
   positions of atoms (default=-1) \param fixed         		 for
   some groups, two different angle dependencies are defined in GRID: one which
   takes some flexibility of groups (rotation/swapping of lone pairs and
   hydrogen) into account and one for strictly fixed conformations if true,
   strictly fixed conformations (default=fixed) \param cutoff
   minimum cutoff distance [A] (default:1.0)
   */
  Hydrophilic(const RDKit::ROMol &mol, int confId = -1, bool fixed = true,
              double cutoff = 1.0);
  ~Hydrophilic() {};

  double operator()(double x, double y, double z, double thres);

 private:
  HBond d_hbondOH, d_hbondO;
};

//! \brief writes the contents of the MIF to a stream
/*
 The Grid \c grd is written in Gaussian Cube format
 A molecule \c mol has to be specified
 */
void writeToCubeStream(const RDGeom::UniformRealValueGrid3D &grd,
                       const RDKit::ROMol &mol, std::ostream &outStrm,
                       int confid = -1);

//! \brief writes the contents of the MIF to a file
/*
 The Grid \c grd is written in Gaussian Cube format
 A molecule \c mol has to be specified
 */
void writeToCubeFile(const RDGeom::UniformRealValueGrid3D &grd,
                     const RDKit::ROMol &mol, const std::string &filename,
                     int confid = -1);

//! \brief reads the contents of the MIF from a stream in Gaussian cube format
/*
 The Grid \c grd is modified according to input values
 A pointer to a molecule is returned with all atoms and its positions but NO
 bond information!
 */
std::unique_ptr<RDKit::RWMol> readFromCubeStream(
    RDGeom::UniformRealValueGrid3D &grd, std::istream &inStrm);

//! \brief reads the contents of the MIF from a file in Gaussian cube format
/*
 The Grid \c grd is modified according to input values
 A pointer to a molecule is returned with all atoms and its positions but NO
 bond information!
 */
std::unique_ptr<RDKit::RWMol> readFromCubeFile(
    RDGeom::UniformRealValueGrid3D &grd, const std::string &filename);

}  // end of namespace RDMIF

#endif /* RDMIF_DESCRIPTORS_H */
