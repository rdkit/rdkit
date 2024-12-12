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
#include <RDGeneral/export.h>
#ifndef RDMIF_DESCRIPTORS_H
#define RDMIF_DESCRIPTORS_H
/*! \file MIFDescriptors.h

  \brief code for generating molecular interaction field (MIF) descriptors

  \b Note that this functionality is experimental and the API and/or results
  amay change in future releases.
*/

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
 (default: 0.5 Angstrom) and a spatial extent of conformer \c confId
 (default: -1) of molecule \c mol plus a defined margin \c margin on each side
 (default: 5.0 Angstrom).

 \param mol		Molecule for which the grid is constructed
 \param confId		Id of Conformer which is used
 \param margin 	distance from molecule to the grid's wall
 \param spacing	grid's spacing

 <b>Returns</b>
 pointer to a grid
 */
RDKIT_MOLINTERACTIONFIELDS_EXPORT
std::unique_ptr<RDGeom::UniformRealValueGrid3D> constructGrid(
    const RDKit::ROMol &mol, int confId = -1, double margin = 5.0,
    double spacing = 0.5);

//! \brief calculates a descriptor at every grid point of MIF
/*!
 any unary function/functor can be used taking a grid point of type
 RDKit::Point3D as parameter

 \param grd			a UniformRealValueGrid3D
 \param functor	any unary function/functor which takes four doubles as
 parameters (3 coordinates, 1 threshold)
 \param thres	  cuts off atom contributions if distance of interaction higher
 than threshold, negative threshold: no threshold, defaults to -1.0
 */
template <typename T>
void calculateDescriptors(RDGeom::UniformRealValueGrid3D &grd, const T &functor,
                          double thres = -1.0) {
  const RDGeom::Point3D &offSet = grd.getOffset();
  auto x = offSet.x;
  auto y = offSet.y;
  auto z = offSet.z;
  auto oX = x;
  auto oY = y;
  auto spacing = grd.getSpacing();
  if (thres < 0) {
    thres = spacing * grd.getSize();
  }
  thres *= thres;  // comparing squared distance

  unsigned int id = 0;
  auto &data = grd.getData();

  for (unsigned int idZ = 0; idZ < grd.getNumZ(); ++idZ) {
    for (unsigned int idY = 0; idY < grd.getNumY(); ++idY) {
      for (unsigned int idX = 0; idX < grd.getNumX(); ++idX) {
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

//! \brief class for calculation of electrostatic interaction (Coulomb energy)
//! between probe and molecule in vacuum (no dielectric)
/*
 d_nAtoms     No of atoms in molecule
 d_softcore   true if softcore interaction is used, else minimum cutoff distance
 is used
 d_probe      charge of probe [e]
 d_absVal     if true, absolute values of interactions are calculated
 d_alpha      softcore interaction parameter [A^2]
 d_cutoff     squared minimum cutoff distance [A^2]
 d_charges    vector of doubles with all partial charges of the atoms in a
 molecule [e] d_pos        vector of doubles with all positions of the atoms in
 a molecule [A] d_prefactor  prefactor taking into account the geometric
 factors,natural constants and conversion of units
 */
class RDKIT_MOLINTERACTIONFIELDS_EXPORT Coulomb {
 public:
  Coulomb() = default;
  ~Coulomb() = default;

  //! \brief constructs Coulomb object from vectors of charges and positions
  /*!
   \param charges     vector of charges [e]
   \param pos         vector of postions [A]
   \param probeCharge charge of probe [e] (default: 1.0)
   \param absVal      if true, negative (favored) values of interactions are
   calculated (default: false)
   \param alpha       softcore interaction parameter [A^2]; if zero,
   a minimum cutoff distance is used (default: 0.0)
   \param cutoff      minimum cutoff distance [A] (default: 1.0)

   */
  Coulomb(const std::vector<double> &charges,
          const std::vector<RDGeom::Point3D> &positions,
          double probeCharge = 1.0, bool absVal = false, double alpha = 0.0,
          double cutoff = 1.0);

  //! \brief constructs Coulomb object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.

   \param mol         molecule object
   \param confId      conformation id which is used to get positions of atoms
   (default=-1)
   \param absVal      if true, absolute values of interactions are
   calculated (default: false)
   \param probeCharge charge of probe [e] (default: 1.0)
   \param prop	     property key for retrieving partial charges
   of atoms (default: "_GasteigerCharge")
   \param alpha       softcore interaction parameter [A^2], if zero, a minimum
   cutoff distance is used (default: 0.0)
   \param cutoff      minimum cutoff distance [A] (default: 1.0)

   */
  Coulomb(const RDKit::ROMol &mol, int confId = -1, double probeCharge = 1.0,
          bool absVal = false, const std::string &prop = "_GasteigerCharge",
          double alpha = 0.0, double cutoff = 1.0);

  //! \brief calculated the electrostatic interaction at point \c pt in the
  //! molecules field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.

   \return electrostatic interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres) const;

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
 distance is used
 d_dielectric: factor of dielectric constants
 d_probe:      charge of probe [e]
 d_absVal:     if true, negative (favored) values of interactions
 are calculated (default: false)
 d_epsilon:    relative permittivity of solvent
 d_xi:         relative permittivity of solute
 d_alpha:      softcore interaction parameter [A^2]
 d_cutoff:     squared minimum cutoff distance [A^2]
 d_charges:    vector of doubles with all partial charges of the atoms
 in a molecule [e]
 d_pos:        vector of Point3Ds with all positions of the atoms
 in a molecule [A]
 d_prefactor:  prefactor taking into account the geometric factors,
 natural constants and conversion of units
 */
class RDKIT_MOLINTERACTIONFIELDS_EXPORT CoulombDielectric {
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
   \param probeCharge charge of probe [e] (default: 1.0)
   \param absVal      if true, negative (favored) values of interactions are
   calculated (default: false)
   \param alpha       softcore interaction parameter [A^2]; if zero,
   a minimum cutoff distance is used (default: 0.0)
   \param cutoff      minimum cutoff distance [A] (default: 1.0)
   \param epsilon     relative permittivity of solvent (default: 80.0)
   \param xi          relative permittivity of solute (default: 4.0)

   */
  CoulombDielectric(const std::vector<double> &charges,
                    const std::vector<RDGeom::Point3D> &positions,
                    double probeCharge = 1.0, bool absVal = false,
                    double alpha = 0.0, double cutoff = 1.0,
                    double epsilon = 80.0, double xi = 4.0);

  //! \brief constructs Coulomb object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.

   \param mol         molecule object
   \param confId      conformation id which is used to get positions of atoms
   (default=-1)
   \param probeCharge charge of probe [e] (default: 1.0)
   \param absVal      if true, negative (favored) values of interactions are
   calculated (default: false)
   \param prop	       property key for retrieving partial charges of atoms
   (default: "_GasteigerCharge")
   \param alpha       softcore interaction parameter [A^2], if zero, a minimum
   cutoff distance is used (default: 0.0)
   \param cutoff      minimum cutoff distance [A] (default: 1.0)
   \param epsilon     relative permittivity of solvent (default: 80.0)
   \param xi          relative permittivity of solute (default: 4.0)
   */
  CoulombDielectric(const RDKit::ROMol &mol, int confId = -1,
                    double probeCharge = 1.0, bool absVal = false,
                    const std::string &prop = "_GasteigerCharge",
                    double alpha = 0.0, double cutoff = 1.0,
                    double epsilon = 80.0, double xi = 4.0);

  //! \brief returns the electrostatic interaction at point \c pt in the
  //! molecule's field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared threshold distance

   \return electrostatic interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres) const;

 private:
  unsigned int d_nAtoms;
  bool d_softcore, d_absVal;
  double d_cutoff, d_probe;
  double d_epsilon, d_xi, d_alpha;
  double d_dielectric;
  std::vector<double> d_charges, d_sp;
  std::vector<double> d_pos;
  // used as a cache
  mutable std::vector<double> d_dists;
};

//! \brief Abstract class for calculation of Van der Waals interaction between
//! probe and molecule at gridpoint \c pt
/*
 * Either the MMFF94 or the UFF VdW term can be calculated with a defined probe
 * atom and molecule.
 * d_cutoff: minimum cutoff distance [A]
 * d_nAtoms: No of atoms in molecule
 * d_R_star_ij: vector of Lennard Jones parameters for every
 * interaction (vdW-radii)
 * d_wellDepth: vector of Lennard Jones parameters for
 * every interaction (well depths)
 * d_pos: vector of positions of atoms in  molecule
 */
class RDKIT_MOLINTERACTIONFIELDS_EXPORT VdWaals {
 public:
  VdWaals() = default;
  VdWaals(const RDKit::ROMol &mol, int confId = -1, double cutoff = 1.0);
  VdWaals(const VdWaals &other) = delete;
  VdWaals &operator=(const VdWaals &other) = delete;
  VdWaals(VdWaals &&other) = default;
  VdWaals &operator=(VdWaals &&other) = default;
  virtual ~VdWaals() = default;

  //! \brief returns the VdW interaction at point \c pt in the molecules field
  //! in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.

   \return vdW interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres) const;

 protected:
  void fillVectors();
  virtual double calcEnergy(double, double, double) const = 0;
  virtual void fillVdwParamVectors(unsigned int atomIdx) = 0;
  double d_cutoff = 1.0;
  unsigned int d_nAtoms = 0;
  std::vector<double> d_pos;
  std::vector<double> d_R_star_ij;
  std::vector<double> d_wellDepth;
  std::unique_ptr<RDKit::ROMol> d_mol;
};

class RDKIT_MOLINTERACTIONFIELDS_EXPORT MMFFVdWaals : public VdWaals {
 public:
  //! \brief constructs VdWaals object which uses MMFF94 from a molecule object
  /*!
  \param mol           molecule object
  \param confId        conformation id which is used to get positions of atoms
  (default: -1)
  \param probeAtomType MMFF94 atom type for the probe atom
  (default: 6, sp3 oxygen)
  \param cutoff        minimum cutoff distance [A] (default: 1.0)
  \param scaling       scaling of VdW parameters to take hydrogen bonds into
  account (default: false)
  */
  MMFFVdWaals(const RDKit::ROMol &mol, int confId = -1,
              unsigned int probeAtomType = 6, bool scaling = false,
              double cutoff = 1.0);
  MMFFVdWaals(const MMFFVdWaals &other) = delete;
  MMFFVdWaals &operator=(const MMFFVdWaals &other) = delete;
  MMFFVdWaals(MMFFVdWaals &&other) = default;
  MMFFVdWaals &operator=(MMFFVdWaals &&other) = default;
  ~MMFFVdWaals() = default;

 private:
  double calcEnergy(double, double, double) const;  // MMFF energy function
  void fillVdwParamVectors(unsigned int atomIdx);
  bool d_scaling;
  std::unique_ptr<RDKit::MMFF::MMFFMolProperties> d_props;
  const ForceFields::MMFF::MMFFVdWCollection *d_mmffVdW;
  const ForceFields::MMFF::MMFFVdW *d_probeParams;
};

class RDKIT_MOLINTERACTIONFIELDS_EXPORT UFFVdWaals : public VdWaals {
 public:
  //! \brief constructs VdWaals object which uses UFF from a molecule object
  /*!
  \param mol           molecule object
  \param confId        conformation id which is used to get positions of atoms
  (default: -1)
  \param probeAtomType UFF atom type for the probe atom (e.g. "O_3", i.e. a
  tetrahedral oxygen)
  \param cutoff        minimum cutoff distance [A] (default: 1.0)
  */
  UFFVdWaals(const RDKit::ROMol &mol, int confId = -1,
             const std::string &probeAtomType = "O_3", double cutoff = 1.0);
  UFFVdWaals(const UFFVdWaals &other) = delete;
  UFFVdWaals &operator=(const UFFVdWaals &other) = delete;
  UFFVdWaals(UFFVdWaals &&other) = default;
  UFFVdWaals &operator=(UFFVdWaals &&other) = default;
  ~UFFVdWaals() = default;

 private:
  double calcEnergy(double, double, double) const;  // UFF energy function
  void fillVdwParamVectors(unsigned int atomIdx);
  const ForceFields::UFF::ParamCollection *d_uffParamColl;
  const ForceFields::UFF::AtomicParams *d_probeParams;
  RDKit::UFF::AtomicParamVect d_params;
};

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
 target donating, the X-H bond direction is saved
 d_plane:        vector of plane vectors in which the lone pairs are
 */
class RDKIT_MOLINTERACTIONFIELDS_EXPORT HBond {
 public:
  HBond() : d_cutoff(1.0), d_probetype(O), d_nInteract(0) {};

  //! \brief constructs HBond object from a molecule object
  /*!
   The molecule \c mol needs to have partial charges set as property of atoms.
   \param mol           molecule object
   \param confId        conformation id which is used to get positions of atoms
   (default=-1)
   \param probeAtomType atom type for the probe atom ("OH", "O", "NH", "N")
   \param fixed         for some groups, two different angle
   dependencies are defined in GRID: one which takes some flexibility of groups
   (rotation/swapping of lone pairs and hydrogen) into account and one for
   strictly fixed conformations if true, strictly fixed conformations
   (default: fixed)
   \param cutoff        minimum cutoff distance [A] (default: 1.0)
   */
  HBond(const RDKit::ROMol &mol, int confId = -1,
        const std::string &probeAtomType = "OH", bool fixed = true,
        double cutoff = 1.0);
  ~HBond() = default;

  //! \brief returns the hydrogen bonding interaction at point \c pt in the
  //! molecules field in [kJ mol^-1]
  /*!
   \param x, y, z     coordinates at which the interaction is calculated
   \param thres       squared max distance until interactions are calc.

   \return hydrogen bonding interaction energy in [kJ mol^-1]
   */
  double operator()(double x, double y, double z, double thres) const;

  unsigned int getNumInteractions() const { return d_nInteract; }

 private:
  boost::uint8_t d_DAprop;

  enum atomtype {
    N,
    O
  };
  double d_cutoff;
  atomtype d_probetype;
  unsigned int d_nInteract;  // number of HBond interactions

  std::vector<atomtype> d_targettypes;
  std::vector<double> d_pos;        // vector of positions of target atoms
  std::vector<double> d_direction;  // vector of sum of lone pair vectors
  std::vector<double> d_lengths;    // vector of lengths of direction vectors
  std::vector<double> d_plane;      // vector of lone pair plane vectors
  // these two are used as a computational cache during operator()
  mutable std::vector<double>
      d_eneContrib;  // energy contributions of all interactions
  mutable std::vector<double>
      d_vectTargetProbe;  // hydrogen bond direction of probe

  // constructing functions
  unsigned int findSpecials(const RDKit::ROMol &mol, int confId, bool fixed,
                            std::vector<unsigned int> &specials);
  unsigned int findAcceptors(const RDKit::ROMol &mol, int confId,
                             const std::vector<unsigned int> &specials);
  unsigned int findDonors(const RDKit::ROMol &mol, int confId,
                          const std::vector<unsigned int> &specials);
  unsigned int findAcceptorsUnfixed(const RDKit::ROMol &mol, int confId,
                                    const std::vector<unsigned int> &specials);
  unsigned int findDonorsUnfixed(const RDKit::ROMol &mol, int confId,
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
class RDKIT_MOLINTERACTIONFIELDS_EXPORT Hydrophilic {
 public:
  //! \brief constructs Hydrophilic object from a molecule object
  /*!
   params:
   \param mol           		 molecule object
   \param confId        		 conformation id which is used to get
   positions of atoms (default: -1)
   \param fixed         		 for some groups, two different angle
   dependencies are defined in GRID: one which takes some flexibility of groups
   (rotation/swapping of lone pairs and hydrogen) into account and one for
   strictly fixed conformations if true, strictly fixed conformations
   (default: fixed)
   \param cutoff  minimum cutoff distance [A] (default: 1.0)
   */
  Hydrophilic(const RDKit::ROMol &mol, int confId = -1, bool fixed = true,
              double cutoff = 1.0);
  ~Hydrophilic() {};

  double operator()(double x, double y, double z, double thres) const;

 private:
  HBond d_hbondOH, d_hbondO;
};

//! \brief writes the contents of the MIF to a stream
/*!
 The Grid \c grd is written in Gaussian Cube format
 A molecule \c mol and a \c confId can optionally be provided
 */
RDKIT_MOLINTERACTIONFIELDS_EXPORT void writeToCubeStream(
    const RDGeom::UniformRealValueGrid3D &grd, std::ostream &outStrm,
    const RDKit::ROMol *mol = nullptr, int confid = -1);

//! \brief writes the contents of the MIF to a file
/*!
 The Grid \c grd is written in Gaussian Cube format
 A molecule \c mol and a \c confId can optionally be provided
 */
RDKIT_MOLINTERACTIONFIELDS_EXPORT void writeToCubeFile(
    const RDGeom::UniformRealValueGrid3D &grd, const std::string &filename,
    const RDKit::ROMol *mol = nullptr, int confid = -1);

//! \brief reads the contents of the MIF from a stream in Gaussian cube format
/*!
 The Grid \c grd is modified according to input values
 If a molecule was associated to the grid on write, a non-null pointer to a
 molecule is returned with atoms and a conformer, but NO bond information.
 If there is no atom information in the cube file, a null pointer is returned.
 */
RDKIT_MOLINTERACTIONFIELDS_EXPORT std::unique_ptr<RDKit::RWMol>
readFromCubeStream(RDGeom::UniformRealValueGrid3D &grd, std::istream &inStrm);

//! \brief reads the contents of the MIF from a file in Gaussian cube format
/*!
 The Grid \c grd is modified according to input values
 If a molecule was associated to the grid on write, a non-null pointer to a
 molecule is returned with atoms and a conformer, but NO bond information.
 If there is no atom information in the cube file, a null pointer is returned.
 */
RDKIT_MOLINTERACTIONFIELDS_EXPORT std::unique_ptr<RDKit::RWMol>
readFromCubeFile(RDGeom::UniformRealValueGrid3D &grd,
                 const std::string &filename);

}  // end of namespace RDMIF

#endif /* RDMIF_DESCRIPTORS_H */
