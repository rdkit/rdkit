//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cmath>
#include <boost/format.hpp>

#include "../RDKitBase.h"
#include "../../RDGeneral/types.h"
#include "../../Geometry/point.h"
#include "StructChecker.h"
#include "Utilites.h"
#include "Stereo.h"

namespace RDKit {
namespace StructureCheck {

static const double PI = 3.14159265359;
static const double ANGLE_EPSILON = (5.0 * PI / 180.);  // 5 degrees
static const double EPS = 0.0000001;  // float type precision ???
                                      /*
                                          // constants for bond definitions
                                      //#define  CIS_TRANS_EITHER  0x03
                                      //#define  CIS_TRANS_SWAPPED 0x08
                                      */
struct npoint_t {
  double x, y, z;
  int number;
};

static double Angle(double x1, double y1, double x2, double y2) {
  // Returns the angle between the two vectors (x1,y1) and (x2,y2) at (0,0).
  double l1, l2;
  double cos_alpha, sin_alpha;
  double result;

  l1 = sqrt(x1 * x1 + y1 * y1);
  l2 = sqrt(x2 * x2 + y2 * y2);
  if (l1 < 0.00001 || l2 < 0.00001) return (0.0);

  cos_alpha = (x1 * x2 + y1 * y2) / (l1 * l2);
  if (cos_alpha > 1.0)  // safeguard against round off errors
    cos_alpha = 1.0;
  else if (cos_alpha < -1.0)
    cos_alpha = -1.0;
  sin_alpha = (x1 * y2 - x2 * y1) / (l1 * l2);

  result = acos(cos_alpha);
  if (sin_alpha < 0.0) result = 2 * PI - result;
  return result;
}

static double Volume(struct npoint_t tetra[4]) {
  // Computes the signed volume of the tetrahedron defined by the four points in
  // tetra[].

  double ax, ay, az, bx, by, bz, cx, cy, cz;

  ax = tetra[1].x - tetra[0].x;
  ay = tetra[1].y - tetra[0].y;
  az = tetra[1].z - tetra[0].z;
  bx = tetra[2].x - tetra[0].x;
  by = tetra[2].y - tetra[0].y;
  bz = tetra[2].z - tetra[0].z;
  cx = tetra[3].x - tetra[0].x;
  cy = tetra[3].y - tetra[0].y;
  cz = tetra[3].z - tetra[0].z;

  return (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
          az * (bx * cy - by * cx));
}

int DubiousStereochemistry(RWMol &mol) {
  /*
   * Checks if there is some ill-defined stereochemistry in the
   * molecule *mp. The function returns a bit set integer which defines
   * the problems encountered.
   */

  int result = 0;
  //   int is_allene, ndb, jatom;

  std::vector<Neighbourhood> neighbour_array;

  SetupNeighbourhood(mol, neighbour_array);

  // look for EITHER bonds
  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond &bond = *mol.getBondWithIdx(i);
    if (RDKit::Bond::UNKNOWN == bond.getBondDir())  //== EITHER
    {
      std::cerr << "bond " << i << ": Either bond found" << std::endl;
      result |= EITHER_BOND_FOUND;
    }
  }
  // look for stereo bonds to non-stereogenic atoms
  for (unsigned i = 0; i < neighbour_array.size(); i++) {
    const Neighbourhood &nbp = neighbour_array[i];
    unsigned nmulti = 0;
    bool is_allene = false;
    for (unsigned j = 0; j < nbp.Bonds.size(); j++) {
      const Bond &bond = *mol.getBondWithIdx(nbp.Bonds[j]);
      if (RDKit::Bond::SINGLE != bond.getBondType()) {
        unsigned ndb = 0;
        nmulti++;
        if (RDKit::Bond::DOUBLE == bond.getBondType()) {
          unsigned jatom = nbp.Atoms[j];
          for (unsigned int jj : neighbour_array[jatom].Bonds)
            if (RDKit::Bond::DOUBLE == mol.getBondWithIdx(jj)->getBondType())
              ndb++;
        }
        if (2 == ndb) is_allene = true;
      }
    }
    unsigned element = mol.getAtomWithIdx(i)->getAtomicNum();
    unsigned n_ligands = (unsigned)nbp.Bonds.size();

    if (!((6 == element &&  // "C"
           n_ligands > 2 && n_ligands <= 4 && nmulti == 0) ||
          (6 == element &&  // "C"
           n_ligands >= 2 && n_ligands <= 3 && nmulti == 1 && is_allene) ||
          (16 == element &&  // "S"
           n_ligands > 2 && n_ligands <= 4) ||
          (7 == element &&  // "N"
           n_ligands > 3 && n_ligands <= 4 && nmulti == 0) ||
          (15 == element &&  // "P"
           n_ligands > 2 && n_ligands <= 4) ||
          ((14 == element &&  // "Si"
            n_ligands > 2 && n_ligands <= 4) &&
           nmulti == 0)))

      for (unsigned j = 0; j < n_ligands; j++) {
        const Bond &bj = *mol.getBondWithIdx(nbp.Bonds[j]);
        if (bj.getBeginAtomIdx() == i &&
            (RDKit::Bond::BEGINWEDGE == bj.getBondDir() ||  // == UP
             RDKit::Bond::BEGINDASH == bj.getBondDir())) {  // == DOWN))
          std::string name;
          mol.getPropIfPresent(common_properties::_Name, name);
          BOOST_LOG(rdWarningLog)
              << boost::format("%10s    atom %3d : %s") % name % i %
                     "stereobond to non-stereogenic atom"
              << std::endl;
          result |= STEREO_BOND_AT_NON_STEREO_ATOM;
        }
      }
  }
  return result;
}

int FixDubious3DMolecule(RWMol &mol) {
  /*
   * Checks if the structure has 3D coordinates and/or flat sp3-carbons with
   * stereo-bonds and
   * converts the designation to 2D, clearing away any Z-component of the
   * coordinates.
   * Real 3D structures without stereo designations go through untouched.
   */

  int result = 0;
  bool non_zero_z = false;
  unsigned nstereo = 0;
  std::vector<Neighbourhood> neighbour_array;
  std::vector<RDGeom::Point3D> atomPoint(
      mol.getNumAtoms());  // X,Y,Z coordinates of each atom

  // X,Y,Z coordinates of each atom
  non_zero_z = getMolAtomPoints(mol, atomPoint);
  // At first check if this is a trivial case i.e. designated '2D'
  // and count the number of stereo bonds
  if (!non_zero_z)  // check Z coordinate of each atom
    return 0;

  nstereo = 0;
  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    if (RDKit::Bond::BEGINWEDGE == bond->getBondDir() ||
        RDKit::Bond::BEGINDASH == bond->getBondDir()
        //???    || RDKit::Bond::EITHERDOUBLE == bond->getBondDir()
    )
      nstereo++;
  }
  if (0 == nstereo) return 0;

  // compute average bond length to use in Volume significance testing
  double length = 0.0;
  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    unsigned a0 = bond->getBeginAtomIdx();
    unsigned a1 = bond->getEndAtomIdx();

    length += (atomPoint[a0].x - atomPoint[a1].x) *
                  (atomPoint[a0].x - atomPoint[a1].x) +
              (atomPoint[a0].y - atomPoint[a1].y) *
                  (atomPoint[a0].y - atomPoint[a1].y) +
              (atomPoint[a0].z - atomPoint[a1].z) *
                  (atomPoint[a0].z - atomPoint[a1].z);
  }
  length = sqrt(length / mol.getNumBonds());

  // check if there is a flat sp3 carbon
  SetupNeighbourhood(mol, neighbour_array);
  int nflat_sp3 = 0;
  for (unsigned i = 0; i < neighbour_array.size(); i++) {
    if (neighbour_array[i].Atoms.size() < 3) continue;
    RDKit::Atom *atom = mol.getAtomWithIdx(i);
    unsigned element = atom->getAtomicNum();
    if (6 != element &&   // "C"
        7 == element &&   // "N"
        15 == element &&  // "P"
        16 == element)    // "S"
      continue;

    unsigned j;
    for (j = 0; j < mol.getNumBonds(); j++) {
      const Bond *bond = mol.getBondWithIdx(j);
      if (RDKit::Bond::BEGINWEDGE == bond->getBondDir() ||
          (RDKit::Bond::BEGINDASH == bond->getBondDir() &&
           i == bond->getBeginAtomIdx()))
        break;
    }
    if (j < mol.getNumBonds()) continue;  // no stereo designation

    double vol = 0.0;
    int stereo_triple;
    const Neighbourhood &nbp = neighbour_array[i];
    unsigned n_ligands = (unsigned)nbp.Bonds.size();
    unsigned i1;
    struct npoint_t tetra[4];
    tetra[0].x = atomPoint[i].x;
    tetra[0].y = atomPoint[i].y;
    tetra[0].z = atomPoint[i].z;
    for (i1 = 0; i1 < n_ligands; i1++)
      if (mol.getBondWithIdx(nbp.Bonds[i1])->getBondType() !=
              RDKit::Bond::SINGLE &&
          16 != element && 15 != element)  // "S" "P"
        break;
    if (i1 >= n_ligands) continue;  // multiple bond found => no sp3 carbon
    stereo_triple = 0;
    for (i1 = 0; i1 < n_ligands; i1++) {
      tetra[1].x = atomPoint[i1].x;
      tetra[1].y = atomPoint[i1].y;
      tetra[1].z = atomPoint[i1].z;
      if (mol.getBondWithIdx(nbp.Bonds[i1])->getBondDir() ==
              RDKit::Bond::BEGINWEDGE ||
          mol.getBondWithIdx(nbp.Bonds[i1])->getBondDir() ==
              RDKit::Bond::BEGINDASH)  // UP DOWN
        stereo_triple |= 1;
      unsigned i2;
      for (i2 = i1 + 1; i2 < n_ligands; i2++) {
        tetra[2].x = atomPoint[i2].x;
        tetra[2].y = atomPoint[i2].y;
        tetra[2].z = atomPoint[i2].z;
        if (mol.getBondWithIdx(nbp.Bonds[i2])->getBondDir() ==
                RDKit::Bond::BEGINWEDGE ||
            mol.getBondWithIdx(nbp.Bonds[i2])->getBondDir() ==
                RDKit::Bond::BEGINDASH)  // UP DOWN
          stereo_triple |= 2;
        unsigned i3;
        for (i3 = i2 + 1; i3 < n_ligands; i3++) {
          tetra[3].x = atomPoint[i3].x;
          tetra[3].y = atomPoint[i3].y;
          tetra[3].z = atomPoint[i3].z;
          if (mol.getBondWithIdx(nbp.Bonds[i3])->getBondDir() ==
                  RDKit::Bond::BEGINWEDGE ||
              mol.getBondWithIdx(nbp.Bonds[i3])->getBondDir() ==
                  RDKit::Bond::BEGINDASH)  // UP DOWN
            stereo_triple |= 4;
          vol = Volume(tetra);
          if (vol < 0.) vol = -vol;
          if (!stereo_triple) continue;
          if (vol < 0.01 * length * length * length) {
            nflat_sp3++;
            break;
          }
          stereo_triple &= ~4;
        }
        stereo_triple &= ~2;
        if (i3 >= n_ligands) break;
      }
      stereo_triple &= ~1;
      if (i2 >= n_ligands) break;
    }
  }

  if (non_zero_z && 0 == mol.getConformer().is3D())  // && mol dim is 2D
    for (unsigned i = 0; i < mol.getNumAtoms(); i++)
      // TODO: ???
      atomPoint[i].z = 0.0;  // set in mol !!!
  result |= ZEROED_Z_COORDINATES;
  // Cleared z-coordinates in 2D MOL file
  if (non_zero_z && nstereo > 0 && nflat_sp3 > 0) {
    RDKit::MolOps::removeStereochemistry(mol);
    result |= CONVERTED_TO_2D;
  }
  return result;
}

void RemoveDubiousStereochemistry(RWMol &mol) {
  // Removes ill-defined stereodescriptors.
  std::vector<Neighbourhood> neighbour_array;
  SetupNeighbourhood(mol, neighbour_array);

  // remove EITHER marks
  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    Bond *bond = mol.getBondWithIdx(i);
    if (RDKit::Bond::UNKNOWN == bond->getBondDir())  //== EITHER
      bond->setBondDir(RDKit::Bond::NONE);
  }
  // remove stereo marks to non-stereogenic atoms
  for (unsigned i = 0; i < neighbour_array.size(); i++) {
    const Neighbourhood &nbp = neighbour_array[i];
    unsigned nmulti = 0;
    for (unsigned j = 0; j < nbp.Atoms.size(); j++) {
      const Bond &bond = *mol.getBondWithIdx(nbp.Bonds[j]);
      if (RDKit::Bond::SINGLE != bond.getBondType()) nmulti++;
    }

    unsigned element = mol.getAtomWithIdx(i)->getAtomicNum();
    unsigned n_ligands = (unsigned)nbp.Bonds.size();
    if (!((6 == element &&  // "C"
           n_ligands > 2 && n_ligands <= 4 && nmulti == 0) ||
          (16 == element &&  // "S"
           n_ligands > 2 && n_ligands <= 4) ||
          (7 == element &&  // "N"
           n_ligands > 3 && n_ligands <= 4 && nmulti == 0) ||
          (15 == element &&  // "P"
           n_ligands > 2 && n_ligands <= 4) ||
          ((14 == element &&  // "Si"
            n_ligands > 2 && n_ligands <= 4) &&
           nmulti == 0))) {
      for (unsigned j = 0; j < n_ligands; j++) {
        Bond &bj = *mol.getBondWithIdx(nbp.Bonds[j]);
        if (bj.getBeginAtomIdx() == i &&
            (RDKit::Bond::BEGINWEDGE == bj.getBondDir()      // == UP
             || RDKit::Bond::BEGINDASH == bj.getBondDir()))  // == DOWN))
          bj.setBondDir(RDKit::Bond::NONE);
      }
    }
  }
}

//----------------------------------------------------------------------
// CheckStereo():
//----------------------------------------------------------------------

struct stereo_bond_t {
  double x, y;                     // relative 2D coordinates
  RDKit::Bond::BondDir direction;  // stereo direction to this atom coord
  int number;                      // atom number of ligand atom
  double angle;                    // angle in radiants rel. to first bond
                                   // in array (counted counter clockwise)
};

static int Atom3Parity(struct stereo_bond_t ligands[3]) {
  // Computes the stereo parity defined by three ligands.
  int a, b;
  int reference;
  struct npoint_t tetrahedron[4], h;
  double angle;
  int maxnum;

  maxnum = ligands[0].number;
  for (unsigned i = 1; i < 3; i++) {
    if (maxnum < ligands[i].number) maxnum = ligands[i].number;
  }

  reference = (-1);
  for (unsigned i = 0; i < 3; i++) {
    if (ligands[i].direction != RDKit::Bond::NONE) {
      if (reference == (-1)) {
        reference = i;
      } else {
        // stereo_error = "three attachments with more than 2 stereobonds";
        std::cerr << "three attachments with more than 2 stereobonds"
                  << std::endl;
        return (ILLEGAL_REPRESENTATION);
      }
    }
  }

  if (reference == (-1)) {
    return (UNDEFINED_PARITY);
  }

  if (reference == 0) {
    a = 1;
    b = 2;
  } else if (reference == 1) {
    a = 0;
    b = 2;
  } else {
    a = 0;
    b = 1;
  }

  angle = Angle(ligands[a].x, ligands[a].y, ligands[b].x, ligands[b].y);

  if (angle < ANGLE_EPSILON || fabs(PI - angle) < ANGLE_EPSILON) {
    // stereo_error = "three attachments: colinearity violation";
    std::cerr << "three attachments colinearity violation" << std::endl;

    return (ILLEGAL_REPRESENTATION);
  }

  tetrahedron[0].x = 0.0;
  tetrahedron[0].y = 0.0;
  tetrahedron[0].z = 0.0;
  tetrahedron[0].number = maxnum + 1;
  for (unsigned i = 0; i < 3; i++) {
    tetrahedron[i + 1].x = ligands[i].x;
    tetrahedron[i + 1].y = ligands[i].y;
    if (ligands[i].direction == RDKit::Bond::BEGINWEDGE)  // UP)
      tetrahedron[i + 1].z = 1.0;
    else if (ligands[i].direction == RDKit::Bond::BEGINDASH)  // DOWN)
      tetrahedron[i + 1].z = -1.0;
    else if (ligands[i].direction == RDKit::Bond::NONE)
      tetrahedron[i + 1].z = 0.0;
    else {
      // stereo_error = "three attachments: illegal bond symbol";
      std::cerr << "three attachments illegal bond symbol" << std::endl;
      return (ILLEGAL_REPRESENTATION);
    }
    tetrahedron[i + 1].number = ligands[i].number;
  }

  for (unsigned i = 1; i < 4; i++)
    for (unsigned j = i; j > 0; j--)
      if (tetrahedron[j].number < tetrahedron[j - 1].number) {
        h = tetrahedron[j];
        tetrahedron[j] = tetrahedron[j - 1];
        tetrahedron[j - 1] = h;
      } else
        break;

  return (Volume(tetrahedron) > 0.0 ? EVEN_PARITY : ODD_PARITY);
}

static int Atom4Parity(struct stereo_bond_t ligands[4]) {
  /*
   * Computes the stereo parity defined by four ligands.
   * Assumes central atom at 0/0/0.
   */
  struct npoint_t tetrahedron[4], h;
  int nup, ndown, nopposite;
  double angle;

  nup = ndown = 0;
  for (unsigned i = 0; i < 4; i++) {
    tetrahedron[i].x = ligands[i].x;
    tetrahedron[i].y = ligands[i].y;
    tetrahedron[i].z = 0.0;
    tetrahedron[i].number = ligands[i].number;
    if (ligands[i].direction == RDKit::Bond::BEGINWEDGE) {  // UP
      nup++;
      tetrahedron[i].z = 1.0;
    } else if (ligands[i].direction == RDKit::Bond::BEGINDASH) {  // DOWN
      ndown++;
      tetrahedron[i].z = (-1.0);
    } else if (ligands[i].direction != RDKit::Bond::NONE) {
      // stereo_error = "illegal bond symbol";
      std::cerr << "illegal bond symbol" << std::endl;

      return (ILLEGAL_REPRESENTATION);
    }
  }

  if (nup == 0 && ndown == 0) return (UNDEFINED_PARITY);

  if (nup > 2 || ndown > 2) {
    // stereo_error = "too many stereobonds";
    std::cerr << "too many stereobonds" << std::endl;
    return (ILLEGAL_REPRESENTATION);
  }

  if (nup + ndown == 1)  // check for 'umbrellas'
  {
    unsigned ij;
    for (ij = 0; ij < 4; ij++)
      if (ligands[ij].direction == RDKit::Bond::BEGINWEDGE ||
          ligands[ij].direction == RDKit::Bond::BEGINDASH)
        break;
    nopposite = 0;
    for (unsigned j = 0; j < 4; j++)
      if (ij == j)
        continue;
      else if (ligands[ij].x * ligands[j].x + ligands[ij].y * ligands[j].y < 0)
        nopposite++;
    if (nopposite > 2) {
      // stereo_error = "UMBRELLA: all non-stereo bonds opposite to single
      // stereo bond";
      std::cerr << "umbrella" << std::endl;
      return (ILLEGAL_REPRESENTATION);
    }
  }

  for (unsigned i = 0; i < 2; i++)
    if ((ligands[i].direction == RDKit::Bond::BEGINWEDGE &&
         ligands[i + 2].direction == RDKit::Bond::BEGINDASH) ||
        (ligands[i].direction == RDKit::Bond::BEGINDASH &&
         ligands[i + 2].direction == RDKit::Bond::BEGINWEDGE)) {
      // stereo_error = "UP/DOWN opposition";
      std::cerr << "up/down" << std::endl;
      return (ILLEGAL_REPRESENTATION);
    }

  for (unsigned i = 0; i < 4; i++)
    if ((ligands[i].direction == RDKit::Bond::BEGINWEDGE &&
         ligands[(i + 1) % 4].direction == RDKit::Bond::BEGINWEDGE)  // UP
        || (ligands[i].direction == RDKit::Bond::BEGINDASH           // DOWN
            && ligands[(i + 1) % 4].direction == RDKit::Bond::BEGINDASH)) {
      // stereo_error = "Adjacent like stereobonds";
      std::cerr << "adjacent like" << std::endl;
      return (ILLEGAL_REPRESENTATION);
    }

  for (unsigned i = 0; i < 4; i++)
    if (ligands[i].direction == RDKit::Bond::NONE &&
        ligands[(i + 1) % 4].direction == RDKit::Bond::NONE &&
        ligands[(i + 2) % 4].direction == RDKit::Bond::NONE) {
      angle = Angle(ligands[i].x - ligands[(i + 1) % 4].x,
                    ligands[i].y - ligands[(i + 1) % 4].y,
                    ligands[(i + 2) % 4].x - ligands[(i + 1) % 4].x,
                    ligands[(i + 2) % 4].y - ligands[(i + 1) % 4].y);
      if (angle < (185 * PI / 180)) {
        // stereo_error = "colinearity or triangle rule violation";
        std::cerr << "colinearity or triangle rule" << std::endl;
        return (ILLEGAL_REPRESENTATION);
      }
    }

  for (unsigned i = 1; i < 4; i++)
    for (unsigned j = i; j > 0; j--)
      if (tetrahedron[j].number < tetrahedron[j - 1].number) {
        h = tetrahedron[j];
        tetrahedron[j] = tetrahedron[j - 1];
        tetrahedron[j - 1] = h;
      } else
        break;

  return (Volume(tetrahedron) > 0.0 ? EVEN_PARITY : ODD_PARITY);
}

int AtomParity(const ROMol &mol, unsigned iatom, const Neighbourhood &nbp) {
  /*
   * Computes the stereo parity of atom number iatom in *mp relative
   * to its numbering. The immediate neighbours are defined by *nbp
   * to speed up processing.
   */
  struct stereo_bond_t stereo_ligands[4], h;
  bool multiple = false, allene = false, stereo = false;

  if (nbp.Atoms.size() < 3 || nbp.Atoms.size() > 4)
    return ILLEGAL_REPRESENTATION;

  std::vector<RDGeom::Point3D> atomPoint(
      mol.getNumAtoms());  // X,Y,Z coordinates of each atom
  getMolAtomPoints(mol, atomPoint);

  for (unsigned i = 0; i < nbp.Bonds.size(); i++) {
    const Bond &bi = *mol.getBondWithIdx(i);
    if (bi.getBondType() != RDKit::Bond::SINGLE) {
      multiple = true;
      // check if the multiple bond is part of an allene
      unsigned jatom = nbp.Atoms[i] + 1;
      unsigned ndb = 0;
      for (unsigned j = 0; j < mol.getNumBonds(); j++) {
        const Bond &bond = *mol.getBondWithIdx(j);
        if (bond.getBeginAtomIdx() == jatom || bond.getEndAtomIdx() == jatom)
          if (bond.getBondType() == RDKit::Bond::DOUBLE) ndb++;
      }
      if (ndb == 2) allene = true;
    }

    stereo_ligands[i].x = atomPoint[i].x - atomPoint[iatom].x;
    stereo_ligands[i].y = atomPoint[i].y - atomPoint[iatom].y;
    stereo_ligands[i].number = nbp.Atoms[i] + 1;
    if (bi.getBeginAtomIdx() == iatom) {
      stereo_ligands[i].direction = bi.getBondDir();
      if (stereo_ligands[i].direction == RDKit::Bond::BEGINWEDGE ||  // UP ||
          stereo_ligands[i].direction == RDKit::Bond::BEGINDASH)     // DOWN
        stereo = true;
    } else
      stereo_ligands[i].direction = RDKit::Bond::NONE;
  }
  unsigned element = mol.getAtomWithIdx(iatom)->getAtomicNum();
  if (multiple && stereo && 15 != element && 16 != element) {  // "P" && "S"
    if (allene)
      return (ALLENE_PARITY);
    else {
      // stereo_error = "AtomParity: Stereobond at unsaturated atom";
      std::cerr << "stereobond at unsaturated atom" << std::endl;
      return (ILLEGAL_REPRESENTATION);
    }
  } else if (multiple && 16 != element)  // "S"
    return (UNDEFINED_PARITY);

  stereo_ligands[0].angle = 0.0; /* comp. angle rel. to first ligand */
  for (unsigned i = 1; i < nbp.Atoms.size(); i++)
    stereo_ligands[i].angle = Angle(stereo_ligands[0].x, stereo_ligands[0].y,
                                    stereo_ligands[i].x, stereo_ligands[i].y);
  for (unsigned i = 2; i < nbp.Atoms.size(); i++) /* sort ligands */
    for (unsigned j = i; j > 1; j--)
      if (stereo_ligands[j].angle < stereo_ligands[j - 1].angle) {
        h = stereo_ligands[j];
        stereo_ligands[j] = stereo_ligands[j - 1];
        stereo_ligands[j - 1] = h;
      } else
        break;

  return (nbp.Atoms.size() == 3 ? Atom3Parity(stereo_ligands)
                                : Atom4Parity(stereo_ligands));
}

bool CheckStereo(const ROMol &mol) {
  /*
   * Checks if all potential stereocenters are either completely undefined
   * or attributed with hashes and wedges according to MDL rules.
   */
  bool result = true;
  int parity;
  bool center_defined = false;
  std::vector<Neighbourhood> neighbour_array;

  SetupNeighbourhood(mol, neighbour_array);

  for (unsigned i = 0; i < neighbour_array.size(); i++) {
    const Neighbourhood &nbp = neighbour_array[i];
    unsigned element = mol.getAtomWithIdx(i)->getAtomicNum();
    unsigned n_ligands = (unsigned)nbp.Bonds.size();
    if ((n_ligands > 2 && n_ligands <= 4) && (6 == element       // "C"
                                              || 16 == element   // "S"
                                              || 7 == element    // "N"
                                              || 8 == element    // "O"
                                              || 15 == element   // "P"
                                              || 14 == element)  // "Si"
    ) {
      parity = AtomParity(mol, i, nbp);

      const Atom *atom = mol.getAtomWithIdx(i);
      if (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)
        center_defined = true;

      if (parity == ILLEGAL_REPRESENTATION) {  // stereo_error
        result = false;
      } else if (parity == EVEN_PARITY || parity == ODD_PARITY ||
                 parity == ALLENE_PARITY) {
        center_defined = true;
      } else {
        for (unsigned j = 0; j < nbp.Bonds.size(); j++) {
          const Bond &bond = *mol.getBondWithIdx(j);
          if (bond.getBeginAtomIdx() == i &&
              (RDKit::Bond::BEGINWEDGE == bond.getBondDir()        // == UP
               || RDKit::Bond::BEGINDASH == bond.getBondDir())) {  // == DOWN))
            // stereobond to non-stereogenic atom
            std::cerr << "stereobond to nonstereogenic" << std::endl;
            result = false;
          }
        }
      }
    }
  }

  unsigned int chiralFlag = 0;
  mol.getPropIfPresent(RDKit::common_properties::_MolFileChiralFlag,
                       chiralFlag);

  if (chiralFlag && !center_defined) {  // no stereocenter defined
    std::cerr << "chiral flag, no stereocenter" << std::endl;
    result = false;
  }
  return result;
}

bool AtomClash(RWMol &mol, double clash_limit) {
  /*
   * Checks if any two atoms in *mp come closer than 10% of the
   * average bond length or if an atom is too close the line
   * between two bonded atoms.
   */
  double bond_square_median, dist, min_dist;
  double rr, rb, bb, h;
  std::vector<RDGeom::Point3D> atomPoint(
      mol.getNumAtoms());  // X,Y,Z coordinates of each atom

  bool twod = true;
  getMolAtomPoints(mol, atomPoint, twod);

  // compute median of square of bond length (quick/dirty)
  if (mol.getNumBonds() == 0) return false;
  std::vector<double> blengths(mol.getNumBonds());
  blengths[0] = 1.0;

  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    unsigned a1 = bond->getBeginAtomIdx();
    unsigned a2 = bond->getEndAtomIdx();
    blengths[i] = (atomPoint[a1].x - atomPoint[a2].x) *
                      (atomPoint[a1].x - atomPoint[a2].x) +
                  (atomPoint[a1].y - atomPoint[a2].y) *
                      (atomPoint[a1].y - atomPoint[a2].y);
  }
  for (unsigned i = 1; i < mol.getNumBonds(); i++) {
    for (int j = rdcast<int>(i) - 1; j >= 0; j--) {
      if (blengths[j] > blengths[j + 1]) {
        h = blengths[j];
        blengths[j] = blengths[j + 1];
        blengths[j + 1] = h;
      } else
        break;
    }
  }
  bond_square_median = blengths[mol.getNumBonds() / 2];

  // Check if two atoms get too close to each other
  min_dist = bond_square_median;
  for (unsigned i = 0; i < mol.getNumAtoms(); i++)
    for (unsigned j = i + 1; j < mol.getNumAtoms(); j++) {
      dist =
          (atomPoint[i].x - atomPoint[j].x) *
              (atomPoint[i].x - atomPoint[j].x) +
          (atomPoint[i].y - atomPoint[j].y) * (atomPoint[i].y - atomPoint[j].y);
      if (dist < clash_limit * clash_limit * bond_square_median) {
        std::cerr << "clash 1" << std::endl;

        return true;
      }
      if (dist < min_dist) min_dist = dist;
    }

  // check if atom lies on top of some bond
  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    unsigned a1 = bond->getBeginAtomIdx();
    unsigned a2 = bond->getEndAtomIdx();
    for (unsigned j = 0; j < mol.getNumAtoms(); j++) {
      if (a1 == j || a2 == j) continue;
      rr = (atomPoint[j].x - atomPoint[a1].x) *
               (atomPoint[j].x - atomPoint[a1].x) +
           (atomPoint[j].y - atomPoint[a1].y) *
               (atomPoint[j].y - atomPoint[a1].y);
      bb = (atomPoint[a2].x - atomPoint[a1].x) *
               (atomPoint[a2].x - atomPoint[a1].x) +
           (atomPoint[a2].y - atomPoint[a1].y) *
               (atomPoint[a2].y - atomPoint[a1].y);
      rb = (atomPoint[j].x - atomPoint[a1].x) *
               (atomPoint[a2].x - atomPoint[a1].x) +
           (atomPoint[j].y - atomPoint[a1].y) *
               (atomPoint[a2].y - atomPoint[a1].y);

      if (0 <= rb &&   // cos alpha > 0
          rb <= bb &&  // projection of r onto b does not exceed b
          (rr * bb - rb * rb) / (bb + EPS) <  // distance from bond < limit
              clash_limit * clash_limit * bond_square_median) {
        BOOST_LOG(rdWarningLog)
            << "atom " << j << " "
            << 100 * sqrt((rr * bb - rb * rb) / (bb * bond_square_median + EPS))
            << "% of average bond length " << a1 << "-" << a2 << std::endl;
        std::cerr << "clash 2" << std::endl;
        return true;
      }
    }
  }
  return false;  // no clash
}
//--------------------------------------------------------------------------

/*
 * Sets the color field of the defined double bonds in *mp to CIS,
 * TRANS, or NONE depending on the ligands with the lowest numbering[].
 * It returns the number of defined double bonds found.
 */
int CisTransPerception(const ROMol &mol,
                       const std::vector<RDGeom::Point3D> &points,
                       const std::vector<unsigned> &numbering,
                       std::vector<unsigned> &bondColor) {
  int result = 0;
  unsigned int maxnum = 0;
  std::vector<Neighbourhood> nba(mol.getNumAtoms());
  SetupNeighbourhood(mol, nba);

  for (unsigned int &i : bondColor) i = 0;
  for (unsigned i = 0; i < nba.size(); i++)  // n_atoms
    if (numbering[i] > maxnum) maxnum = numbering[i];

  for (unsigned i = 0; i < bondColor.size(); i++)
    if (RDKit::Bond::DOUBLE == mol.getBondWithIdx(i)->getBondType()
        // FIX:         && mol.getBondWithIdx(i)->getBondDir() !=
        // RDKit::Bond::ENDDOWNRIGHT
        // FIX:         && mol.getBondWithIdx(i)->getBondDir() !=
        // RDKit::Bond::ENDUPRIGHT
        && mol.getBondWithIdx(i)->getBondDir() !=
               RDKit::Bond::EITHERDOUBLE) {  // != CIS_TRANS_EITHER
      unsigned j1 = mol.getBondWithIdx(i)->getBeginAtomIdx();
      unsigned j2 = mol.getBondWithIdx(i)->getEndAtomIdx();
      if (6 != mol.getAtomWithIdx(j1)->getAtomicNum())  // C
        continue;
      if (16 != mol.getAtomWithIdx(j1)->getAtomicNum())  // N
        continue;
      if (6 != mol.getAtomWithIdx(j2)->getAtomicNum())  // C
        continue;
      if (16 != mol.getAtomWithIdx(j2)->getAtomicNum())  // N
        continue;
      // n_ligands :
      if (nba[j1].Atoms.size() <= 1 ||  // no subst.
          nba[j2].Atoms.size() <= 1)
        continue;
      if (nba[j1].Atoms.size() > 3 ||  // valence error in mol
          nba[j2].Atoms.size() > 3)
        continue;

      bool equal = false;  // find lowest numbered neighbour of j1
      unsigned at1 = 0;
      for (unsigned k = 0, nmin = maxnum; k < nba[j1].Atoms.size(); k++)
        if (nba[j1].Atoms[k] != j2)  // no loop back
        {
          if (numbering[nba[j1].Atoms[k]] <
              nmin) {  // numbering[nba[j1].atoms[k]]
            at1 = nba[j1].Atoms[k];
            nmin = numbering[at1];
          } else if (numbering[nba[j1].Atoms[k]] == nmin)
            equal = true;
        }
      if (equal)   // identical substituents
        continue;  // no stereochemistry

      equal = false;  // find lowest numbered neighbour of j1
      unsigned at2 = 0;
      for (unsigned k = 0, nmin = maxnum; k < nba[j2].Atoms.size(); k++)  //
        if (nba[j2].Atoms[k] != j1) {  // no loop back
          if (numbering[nba[j2].Atoms[k]] < nmin) {
            at2 = nba[j2].Atoms[k];
            nmin = numbering[at2];
          } else if (numbering[nba[j2].Atoms[k]] == nmin)
            equal = true;
        }
      if (equal)   // identical substituents
        continue;  // no stereochemistry

      // Now, bp points to a double bond, at1 and at2 are the
      // indices (not numbers) of the atoms with lowest numbering
      // attached to each end of the bond, and the bond is
      // guaranteed to be either CIS or TRANS

      double x21, y21;
      double x23, y23;
      double x32, y32;
      double x34, y34;
      double sign;
      x21 = points[at1].x - points[j1].x;
      y21 = points[at1].y - points[j1].y;
      x23 = points[j2].x - points[j1].x;
      y23 = points[j2].y - points[j1].y;
      x32 = (-x23);
      y32 = (-y23);
      x34 = points[at2].x - points[j2].x;
      y34 = points[at2].y - points[j2].y;
      sign = (x21 * y23 - x23 * y21) * (x32 * y34 - x34 * y32);
      if (fabs(sign) < 0.001) continue;
      result++;
      bondColor[i] = (sign > 0.0 ? CIS : TRANS);
    }
  return (result);
}

}  // namespace StructureCheck
}  // namespace RDKit
