/*=================================================================*/
/* Copyright (C)  2024  Greg Landrum and other RDKit contributors  */
/* Contributed by NextMove Software, Cambridge, UK.                */
/*                                                                 */
/*                                                                 */
/* @@ All Rights Reserved @@                                       */
/* The contents are covered by the terms of the                    */
/* BSD license, which is included in the file                      */
/* license.txt.                                                    */
/*=================================================================*/
#define _USE_MATH_DEFINES

#include <iostream>
#include <limits>
#include <string>
#include <list>
#include <cmath>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "DCLV.h"
#include "DCLV_dots.h"

using RDGeom::Point3D;

namespace RDKit {
namespace Descriptors {

constexpr int VOXORDER = 16;
int recordCache = -1;

static bool checkExcludedAtoms(const Atom *atm, bool includeLigand) {
  // helper to check whether atom should be included

  const auto *info = atm->getMonomerInfo();
  if (info) {
    const std::string &resName =
        static_cast<const AtomPDBResidueInfo *>(info)->getResidueName();

    switch (resName[0]) {
      case 'D':
        if (resName == "DOD" || resName == "D20"){
          return true;
        }
        break;
      case 'H':
        if (resName == "HOH" || resName == "H20") {
          return true;
        }
        break;
      case 'S':
        if (resName == "SOL" || resName == "SO4" || resName == "SUL")
          return true;
        break;
      case 'W':
        if (resName == "WAT")
          return true;
        break;
      case 'T':
        if (resName == "TIP")
          return true;
        break;
      case 'P':
        if (resName == "P04")
          return true;
        break;
    }

    if (!includeLigand &&
        static_cast<const AtomPDBResidueInfo *>(info)->getIsHeteroAtom()) {
      return true;
    } 
  }

  return false;
}

static bool includeAsPolar(const Atom *atm, const ROMol &mol,
                         bool includeSandP, bool includeHs) {
  // using Peter Ertl definition, polar atoms = O, N, P, S and attached Hs

  switch (atm->getAtomicNum()) {
    // nitrogen and oxygen
    case 7:
    case 8: {
      return true;
    }

    // sulphur and phosphorous
    case 15:
    case 16: {
      return includeSandP;
    }

    // hydrogen
    case 1: {
      if (!includeHs) {
        return false;
      } else {
        for (const auto nbr : mol.atomNeighbors(atm)) {
          if (includeAsPolar(nbr, mol, includeSandP, includeHs)) {
            return true;
          }
        }
      }
      return false;
    }
    // everything else
    default: {
      return false;
    }
  }
}

struct State {
  // upper and lower grid bounds
  double voxX = 0.0;
  double voxY = 0.0;
  double voxZ = 0.0;
  double voxU = 0.0;
  double voxV = 0.0;
  double voxW = 0.0;
  
  std::vector<unsigned int> grid[VOXORDER][VOXORDER][VOXORDER];

  void createVoxelGrid(const Point3D &minXYZ, const Point3D &maxXYZ,
                       const std::vector<Point3D> &positions,
                       const std::vector<double> &radii_) {
    voxX = VOXORDER / ((maxXYZ.x - minXYZ.x) + 0.1);
    voxU = minXYZ.x;
    voxY = VOXORDER / ((maxXYZ.y - minXYZ.y) + 0.1);
    voxV = minXYZ.y;
    voxZ = VOXORDER / ((maxXYZ.z - minXYZ.z) + 0.1);
    voxW = minXYZ.z;

    unsigned int pcount = (unsigned int)positions.size();
    for (unsigned int i = 0; i < pcount; i++) {
      if (radii_[i] != 0.0) {
        // get grid positions and add to list
        const Point3D &pos = positions[i];
        const unsigned int x = (unsigned int)(voxX * (pos.x - voxU));
        const unsigned int y = (unsigned int)(voxY * (pos.y - voxV));
        const unsigned int z = (unsigned int)(voxZ * (pos.z - voxW));
        grid[x][y][z].push_back(i);
      }
    }
  }

  std::vector<std::vector<unsigned int>> findNeighbours(
      const ROMol &mol, const std::vector<Point3D> &positions,
      const std::vector<double> &radii_, double probeRad, double maxRadius) {
    std::vector<std::vector<unsigned int>> nbrs;
    std::vector<unsigned int> atomNeighbours;

    for (const auto atom : mol.atoms()) {
      atomNeighbours.clear();
      const auto atm_idx = atom->getIdx();
      if (radii_[atm_idx] == 0.0) {
        continue;
      }
      const Point3D &pos = positions[atm_idx];
      const double range = radii_[atm_idx] + probeRad + probeRad;
      const double maxDist = range + maxRadius;

      int lx = voxX * (pos.x - maxDist - voxU);
      if (lx < 0) {
        lx = 0;
      }

      int ly = voxY * (pos.y - maxDist - voxV);
      if (ly < 0) {
        ly = 0;
      }
      int lz = voxZ * (pos.z - maxDist - voxW);
      if (lz < 0) {
        lz = 0;
      }

      int ux = voxX * (pos.x + maxDist - voxU);
      if (ux >= VOXORDER) {
        ux = VOXORDER - 1;
      }
      int uy = voxY * (pos.y + maxDist - voxV);
      if (uy >= VOXORDER) {
        uy = VOXORDER - 1;
      }
      int uz = voxZ * (pos.z + maxDist - voxW);
      if (uz >= VOXORDER) {
        uz = VOXORDER - 1;
      }

      for (auto x = lx; x <= ux; x++) {
        for (auto y = ly; y <= uy; y++) {
          for (auto z = lz; z <= uz; z++) {
            for (auto idx : grid[x][y][z]) {
              if (idx != atm_idx) {
                auto dist = range + radii_[idx];
                if ((pos - positions[idx]).lengthSq() < (dist * dist)){
                  atomNeighbours.push_back(idx);
                }
              }
            }                           
          }
        }
      }
      nbrs.push_back(atomNeighbours);
    }
    return nbrs;
  }
};

// constructor definition
DoubleCubicLatticeVolume::DoubleCubicLatticeVolume(const ROMol &mol,
                                                   std::vector<double> radii,
                                                   bool isProtein,
                                                   bool includeLigand,
                                                   double probeRadius,
                                                   int confId)
    : mol(mol),
      radii_(std::move(radii)),
      isProtein(isProtein),
      includeLigand(includeLigand),
      probeRadius(probeRadius),
      confId(confId) {
  //! Class for calculation of the Shrake and Rupley surface area and volume
  //! using the Double Cubic Lattice Method.
  //!
  //! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and
  //! Michael Scharf, "The Double Cubic Lattice Method: Efficient Approaches
  //! to Numerical Integration of Surface Area and Volume and to Dot Surface
  //! Contouring of Molecular Assemblies", Journal of Computational Chemistry,
  //! Vol. 16, No. 3, pp. 273-284, 1995.

  /*!

    \param mol: input molecule or protein
    \param radii: radii for atoms of input mol
    \param isProtein: flag to calculate burried surface area of a protein ligand
    complex [default=false, free ligand]
    \param includeLigand: flag to trigger
    inclusion of bound ligand in surface area and volume calculations where
    molecule is a protein [default=true]
    \param probeRadius: radius of the
    sphere representing the probe solvent atom
    \param confId: conformer ID to consider [default=-1]
    \return class
    object
  */

  positions = mol.getConformer(confId).getPositions();
  maxRadius = *std::max_element(radii_.begin(), radii_.end());

  // total x,y,z centres
  Point3D cXYZ;
  cXYZ.x = cXYZ.y = cXYZ.z = 0.0;

  // total min max coords
  Point3D minXYZ;
  Point3D maxXYZ;

  minXYZ.x = minXYZ.y = minXYZ.z = std::numeric_limits<double>::infinity();
  maxXYZ.x = maxXYZ.y = maxXYZ.z = -std::numeric_limits<double>::infinity();

  State s = State();

  bool init = false;

  for (const auto atom : mol.atoms()) {
    const unsigned int atomIdx = atom->getIdx();
    numAtoms++;

    if (isProtein) {
      if (checkExcludedAtoms(atom, includeLigand)) {
        radii_[atomIdx] = 0.0;
        numAtoms--;
      }
    }

    if (radii_[atomIdx] != 0.0) {
      const Point3D position = positions[atomIdx];
      // get sum over centres
      cXYZ.x += position.x;
      cXYZ.y += position.y;
      cXYZ.z += position.z;

      // loop over atoms to find min + max x, y + z
      if (init) {
        if (position.x > maxXYZ.x) {
          maxXYZ.x = position.x;
        } else if (position.x < minXYZ.x) {
          minXYZ.x = position.x;
        }

        if (position.y > maxXYZ.y) {
          maxXYZ.y = position.y;
        } else if (position.y < minXYZ.y) {
          minXYZ.y = position.y;
        }

        if (position.z > maxXYZ.z) {
          maxXYZ.z = position.z;
        } else if (position.z < minXYZ.z) {
          minXYZ.z = position.z;
        }
      } else {
        maxXYZ.x = minXYZ.x = position.x;
        maxXYZ.y = minXYZ.y = position.y;
        maxXYZ.z = minXYZ.z = position.z;
        init = true;
      }
    }
  }

  // determineCentreOfGravity
  centreOfGravity.x = cXYZ.x / numAtoms;
  centreOfGravity.y = cXYZ.y / numAtoms;
  centreOfGravity.z = cXYZ.z / numAtoms;

  if (init) {
    s.createVoxelGrid(minXYZ, maxXYZ, positions, radii_);
    neighbours =
        s.findNeighbours(mol, positions, radii_, probeRadius, maxRadius);
  }
}

bool DoubleCubicLatticeVolume::testPoint(
    const Point3D &vect, double solvrad,
    const std::vector<unsigned int> &nbrs) {
  if (recordCache != -1) {
    const Point3D &pos = positions[recordCache];
    const double dist = radii_[recordCache] + solvrad;
    if ((pos - vect).lengthSq() < (dist * dist)) {
      return false;
    }
    recordCache = -1;
  }

  for (unsigned int i : nbrs) {
    const Point3D &pos = positions[i];
    const double dist = radii_[i] + solvrad;
    if ((pos - vect).lengthSq() < (dist * dist)) {
      recordCache = i;
      return false;
    }
  }
  return true;
}

double DoubleCubicLatticeVolume::getAtomSurfaceArea(unsigned int atomIdx) {

  // surface area for single atom

  const auto rad = radii_[atomIdx];
  if (rad == 0.0) {
    return 0.0;  // don't include if radius = 0, masked atom
  }

  const Point3D &pos = positions[atomIdx];
  const std::vector<unsigned int> &nbr = neighbours[atomIdx];

  double atomSurfaceArea = 0.0;
  const double factor = rad + probeRadius;

  recordCache = -1;

  // standard dots has fixed NUMDOTS entries
  // using precomputed dots in DCLV_dots.h
  for (int i = 0; i < NUMDOTS; i++) {
    const Point3D &dots =
        Point3D(standardDots[i][0], standardDots[i][1], standardDots[i][2]);
    const auto vect = pos + dots * factor;
    surfacePoints[atomIdx].push_back(vect);

    if (testPoint(vect, probeRadius, nbr)) {
      atomSurfaceArea += dotArea;
    }
  }
  atomSurfaceArea *= ((4.0 * M_PI) * (factor * factor)) / standardArea;

  return atomSurfaceArea;
}

double DoubleCubicLatticeVolume::getSurfaceArea() {
  // check if already calculated in compactness calc
  if (surfaceArea != 0.0) {
    return surfaceArea;
  }

  for (const auto atom : mol.atoms()) {
    const auto atomIdx = atom->getIdx();
    if (radii_[atomIdx] != 0.0) {
      surfaceArea += getAtomSurfaceArea(atomIdx);
    }
  }
  return surfaceArea;
}

double DoubleCubicLatticeVolume::getPartialSurfaceArea(
    const boost::dynamic_bitset<> &incAtoms) {
  // function to get the surface area of a set of atoms
  // input is a set of indexes to include

  double area = 0.0;
  for (const auto atom : mol.atoms()) {
    const auto atomIdx = atom->getIdx();
    bool incAtom = incAtoms[atomIdx];
    if (incAtom && radii_[atomIdx] != 0.0) {
      area += getAtomSurfaceArea(atomIdx);
    }
  }
  return area;
}

double DoubleCubicLatticeVolume::getPolarSurfaceArea(bool includeSandP, bool includeHs) {
  boost::dynamic_bitset<> polarAtoms(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    if (includeAsPolar(atom, mol, includeSandP, includeHs)) {
      polarAtoms.set(atom->getIdx());
    }
  }
  return getPartialSurfaceArea(polarAtoms);
}

std::map<unsigned int, std::vector<RDGeom::Point3D>> &DoubleCubicLatticeVolume::getSurfacePoints() {
  if (!surfacePoints.empty()) {
    return surfacePoints;
  }

  for (const auto atom : mol.atoms()) {
    const auto atomIdx = atom->getIdx();
    if (radii_[atomIdx] != 0.0) {
      const Point3D &pos = positions[atomIdx];
      const double factor = radii_[atomIdx] + probeRadius;

      // standard dots has fixed NUMDOTS entries
      // using precomputed dots in DCLV_dots.h
      for (int i = 0; i < NUMDOTS; i++) {
        const Point3D &dots =
            Point3D(standardDots[i][0], standardDots[i][1], standardDots[i][2]);
        const auto vect = pos + dots * factor;
        surfacePoints[atomIdx].push_back(vect);
      }
    }
  }
  return surfacePoints;
}

double DoubleCubicLatticeVolume::getAtomVolume(unsigned int atomIdx,
                                               double solventRadius) {
  // get volume for single atom

  const double rad = radii_[atomIdx];
  if (rad == 0.0) {
    return 0.0;  // don't include if radius = 0, masked atom
  }

  const auto &pos = positions[atomIdx];
  const auto &nbr = neighbours[atomIdx];

  const auto factor = rad + solventRadius;

  recordCache = -1;

  RDGeom::Point3D p;
  unsigned int count = 0;
  // using precomputed dots in DCLV_dots.h
  for (const auto &dot : standardDots) {
    auto vect = pos + dot * factor;

    if (testPoint(vect, solventRadius, nbr)) {
      p += dot;
      count++;
    }
  }

  /* Calculate Volume with Gauss-Ostrogradskii Theorem */
  auto atomvolume = (pos - centreOfGravity).dotProduct(p);
  atomvolume = factor * factor * (atomvolume + factor * count);
  atomvolume *= ((4.0 / 3) * M_PI) / NUMDOTS;  // 320 standard dots
  return atomvolume;
}

double DoubleCubicLatticeVolume::getVolume() {
  // function to return volume enclosed by probe sphere

  // check if already calculated in compactness or packing density
  if (totalVolume != 0.0) {
    return totalVolume;
  }

  for (const auto atom : mol.atoms()) {
    const unsigned int atomIdx = atom->getIdx();
    if (radii_[atomIdx] != 0.0) {
      totalVolume += getAtomVolume(atomIdx, probeRadius);
    }
  }

  return totalVolume;
}

double DoubleCubicLatticeVolume::getVDWVolume() {
  // check if already calculated in packing density
  if (vdwVolume != 0.0) {
    return vdwVolume;
  }

  // function to return VDW volume (probe radius == 0)
  for (const auto atom : mol.atoms()) {
    const unsigned int atomIdx = atom->getIdx();
    if (radii_[atomIdx] != 0.0) {
      vdwVolume += getAtomVolume(atomIdx, 0.0);
    }
  }

  return vdwVolume;
}

double DoubleCubicLatticeVolume::getPartialVolume(
    const boost::dynamic_bitset<> &incAtoms) {
  // function to get the surface area of a set of atoms
  // input is a set of indexes to include

  double vol = 0.0;
  for (const auto atom : mol.atoms()) {
    const auto atomIdx = atom->getIdx();
    bool incAtom = incAtoms[atomIdx];
    if (incAtom && radii_[atomIdx] != 0.0) {
      vol += getAtomVolume(atomIdx, probeRadius);
    }
  }
  return vol;
}

double DoubleCubicLatticeVolume::getPolarVolume(bool includeSandP, bool includeHs) {
  boost::dynamic_bitset<> polarAtoms(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    if (includeAsPolar(atom, mol, includeSandP, includeHs)) {
      polarAtoms.set(atom->getIdx());
    }
  }
  return getPartialVolume(polarAtoms);
}

double DoubleCubicLatticeVolume::getCompactness() {
  // check both surface area and volume already computed
  if (surfaceArea == 0.0) {
    getSurfaceArea();
  }

  if (totalVolume == 0.0) {
    getVolume();
  }

  return surfaceArea / cbrt(36.0 * M_PI * totalVolume * totalVolume);
}

double DoubleCubicLatticeVolume::getPackingDensity() {
  if (vdwVolume == 0.0) {
    getVDWVolume();
  }

  if (totalVolume == 0.0) {
    getVolume();
  }

  return vdwVolume / totalVolume;
}

}  // namespace Descriptors
}  // namespace RDKit