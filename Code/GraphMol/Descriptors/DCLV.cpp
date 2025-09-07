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

static bool checkExcludedAtoms(const Atom *atm, const bool &includeLigand) {
  // helper to check whether atom should be included
  // radius = 0 if atom is to be excluded

  // check if solvent atom
  const AtomMonomerInfo *info = atm->getMonomerInfo();
  
  if (info){
    const std::string &resName = ((AtomPDBResidueInfo *)info)->getResidueName();
    
    switch (resName[0]) {
      case 'D':
        return resName == "DOD" || resName == "D20";
      case 'H':
        return resName == "HOH" || resName == "H20";
      case 'S':
        return resName == "SOL" || resName == "SO4" || resName == "SUL";
      case 'W':
        return resName == "WAT";
      case 'T':
        return resName == "TIP";
      case 'P':
        return resName == "P04";
      default:
        return false;
    }

    if (!includeLigand && ((AtomPDBResidueInfo *)info)->getIsHeteroAtom()){
      return true;
    }
  }
    
  return false;
}

static bool includeInPSA(const Atom *atm, const ROMol &mol, const bool &includeSandP){
  // using Peter Ertl definition, polar atoms = O, N, P, S and attached Hs
  
  switch(atm->getAtomicNum()){
    // nitrogen and oxygen
    case 7:
    case 8:
    {
      return true;
    }

    // sulphur and phosphorous
    case 15:
    case 16: {
      return includeSandP;
    }
  
    // hydrogen
    case 1: {
      bool isPSAH = false;
      for(const auto &nbri: make_iterator_range(mol.getAtomBonds(atm))) {
        const Bond *bond = (mol)[nbri];
        const Atom *nbr = bond->getOtherAtom(atm);
        if (includeInPSA(nbr, mol, includeSandP)){
          isPSAH = true;
        }
      }
      return isPSAH;
    }
    // everything else
    default: {
      return false;
    }  
  }
}

static bool within(const Point3D &pos1, const Point3D &pos2, const double &dist){
  
  const double dx = pos1.x - pos2.x;
  const double dy = pos1.y - pos2.y;
  const double dz = pos1.z - pos2.z;
  
  return ((dx * dx + dy * dy + dz * dz) < (dist * dist));
}

struct State {
  
  double voxX = 0.0;
  double voxY = 0.0;
  double voxZ = 0.0;
  double voxU = 0.0;
  double voxV = 0.0;
  double voxW = 0.0;

  std::vector<unsigned int> grid[VOXORDER][VOXORDER][VOXORDER];

  void createVoxelGrid(const Point3D &minXYZ, const Point3D &maxXYZ, const std::vector<Point3D> &positions, const std::vector<double> &radii) {

    voxX = VOXORDER / ((maxXYZ.x - minXYZ.x) + 0.1);
    voxU = minXYZ.x;
    voxY = VOXORDER / ((maxXYZ.y - minXYZ.y) + 0.1);
    voxV = minXYZ.y;
    voxZ = VOXORDER / ((maxXYZ.z - minXYZ.z) + 0.1);
    voxW = minXYZ.z;
    
    unsigned int pcount = (unsigned int)positions.size();
    for (unsigned int i = 0; i < pcount; i++) {
      if (radii[i] != 0.0) {
        // get grid positions and add to list
        const Point3D &pos = positions[i];
        const int x = voxX * (pos.x - voxU);
        const int y = voxY * (pos.y - voxV);
        const int z = voxZ * (pos.z - voxW);
        grid[x][y][z].push_back(i);
      }
    }
  }
            
  std::vector<std::vector<unsigned int>> findNeighbours(const ROMol &mol, std::vector<Point3D> &positions, std::vector<double> &radii, double &probeRad){
    
    std::vector<std::vector<unsigned int>> nbrs;
    double maxRadius = 1.87;

    for (const auto atom : mol.atoms()) {
      std::vector<unsigned int> atomNeighbours;
      const unsigned int atm_idx = atom->getIdx();
      if (radii[atm_idx] != 0.0){
        const Point3D pos = positions[atm_idx];
        const double range = radii[atm_idx] + probeRad + probeRad;
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

        for (int x = lx; x <= ux; x++) {
          for (int y = ly; y <= uy; y++) {
            for (int z = lz; z <= uz; z++) {
              std::vector<unsigned int> tmp = grid[x][y][z];
              for (unsigned int idx : tmp) {
                if (idx != atm_idx) {
                  double dist = range + radii[idx];
                  if (within(pos, positions[idx], dist)){
                    atomNeighbours.push_back(idx);
                  }
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
                                                   bool isProtein,
                                                   bool includeLigand,
                                                   double probeRadius,
                                                   int confId) : 
                                                   mol(mol), 
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
  const PeriodicTable *tbl = PeriodicTable::getTable();

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
    const unsigned int aidx = atom->getIdx();
    double rad = tbl->getRvdw(atom->getAtomicNum());

    numAtoms++; 
    
    if (isProtein){
      if (checkExcludedAtoms(atom, includeLigand)){
        rad = 0.0;
        numAtoms--; 
      } 
    }

    if (rad != 0.0){
      const Point3D position = positions[aidx];
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
    radii.push_back(rad); // get radii
  }
  
  // determineCentreOfGravity
  centreOfGravity.x = cXYZ.x / numAtoms;
  centreOfGravity.y = cXYZ.y / numAtoms;
  centreOfGravity.z = cXYZ.z / numAtoms;

  if (init){
    s.createVoxelGrid(minXYZ, maxXYZ, positions, radii);
    neighbours = s.findNeighbours(mol, positions, radii, probeRadius);
  }
}

bool DoubleCubicLatticeVolume::testPoint(double *vect, const double &solvrad, const std::vector<unsigned int> &nbrs){
  
  if (recordCache != -1) {
    const Point3D &pos = positions[recordCache];
    const double dist = radii[recordCache] + solvrad;
    const double dx = pos.x - vect[0];
    const double dy = pos.y - vect[1];
    const double dz = pos.z - vect[2];
    if ((dx * dx + dy * dy + dz * dz) < (dist * dist)) {
      return false;
    }
    recordCache = -1;
  }

  for (unsigned int i : nbrs) {
    const Point3D &pos = positions[i];
    const double dist = radii[i] + solvrad;
    const double dx = pos.x - vect[0];
    const double dy = pos.y - vect[1];
    const double dz = pos.z - vect[2];
    if ((dx * dx + dy * dy + dz * dz) < (dist * dist)) {
      recordCache = i;
      return false;
    }
  }
  return true;
}

double DoubleCubicLatticeVolume::getAtomSurfaceArea(const unsigned int &atom_idx) {
  // surface area for single atom
  
  const double rad = radii[atom_idx];
  if (rad == 0.0){
    return -1.0; // don't include if radius = 0, masked atom
  }

  const Point3D &pos = positions[atom_idx];
  const std::vector<unsigned int> &nbr = neighbours[atom_idx];
  double vect[3];

  double atomSurfaceArea = 0.0;
  const double factor = rad + probeRadius;

  recordCache = -1;

  // standard dots has fixed 320 entries
  // using precomputed dots in DCLV_dots.h
  for (int i = 0; i < 320; i++) {
    
    const Point3D &dots = Point3D(standardDots[i][0], standardDots[i][1], standardDots[i][2]);
    vect[0] = pos.x + factor * dots.x;
    vect[1] = pos.y + factor * dots.y;
    vect[2] = pos.z + factor * dots.z;

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
    const unsigned int atom_idx = atom->getIdx();
    if (radii[atom_idx] != 0.0){ 
      surfaceArea += getAtomSurfaceArea(atom_idx);
    }
  }
  return surfaceArea;
}

double DoubleCubicLatticeVolume::getPolarSurfaceArea(const bool &includeSandP) {
  // TODO check this
  for (const auto atom : mol.atoms()) {
    const unsigned int atom_idx = atom->getIdx();
    bool incAtom = includeInPSA(atom, mol, includeSandP);
    if (radii[atom_idx] != 0.0 && incAtom){ 
      polarSurfaceArea += getAtomSurfaceArea(atom_idx);
    }
  }
  return polarSurfaceArea;
}

double DoubleCubicLatticeVolume::getAtomVolume(const unsigned int &atom_idx, const double &solventRadius) {
  // get volume for single atom
  
  const double rad = radii[atom_idx];
  if (rad == 0.0){
    return -1.0; // don't include if radius = 0, masked atom
  }

  const Point3D &pos = positions[atom_idx];
  const std::vector<unsigned int> &nbr = neighbours[atom_idx];
  double vect[3];

  const double factor = rad + solventRadius;
  unsigned int count = 0;
  double px, py, pz;
  px = py = pz = 0.0;

  recordCache = -1;

  // using precomputed dots in DCLV_dots.h
  for (int i = 0; i < 320; i++) {
    
    const Point3D &dots = Point3D(standardDots[i][0], standardDots[i][1], standardDots[i][2]);
    vect[0] = pos.x + factor * dots.x;
    vect[1] = pos.y + factor * dots.y;
    vect[2] = pos.z + factor * dots.z;

    if (testPoint(vect, solventRadius, nbr)) {
      px += dots.x;
      py += dots.y;
      pz += dots.z;
      count++;
    }
  }

  /* Calculate Volume with Gauss-Ostrogradskii Theorem */
  double atomvolume = (pos.x - centreOfGravity.x) * px + 
                      (pos.y - centreOfGravity.y) * py +
                      (pos.z - centreOfGravity.z) * pz;
  atomvolume = factor * factor * (atomvolume + factor * count);

  atomvolume *= ((4.0 / 3) * M_PI) / 320; // 320 standard dots
  return atomvolume; 
}

double DoubleCubicLatticeVolume::getVolume () {
  
  // function to return volume enclosed by probe sphere

  // check if already calculated in compactness or packing density
  if (totalVolume != 0.0){
    return totalVolume;
  }

  for (const auto atom : mol.atoms()) {
    const unsigned int atom_idx = atom->getIdx();
    if (radii[atom_idx] != 0.0){
      totalVolume += getAtomVolume(atom_idx, probeRadius);
    }
  }

  return totalVolume;
}

double DoubleCubicLatticeVolume::getVDWVolume () {

  // check if already calculated in packing density
  if (vdwVolume != 0.0){
    return vdwVolume;
  }

  // function to return VDW volume (probe radius == 0)
  for (const auto atom : mol.atoms()) {
    const unsigned int atom_idx = atom->getIdx();
    if (radii[atom_idx] != 0.0){
      vdwVolume += getAtomVolume(atom_idx, 0.0);
    }
  }
  
  return vdwVolume;
}

double DoubleCubicLatticeVolume::getCompactness() {
  // check both surface area and volume already computed
  if (surfaceArea == 0.0){
    getSurfaceArea();
  }

  if (totalVolume == 0.0){
    getVolume();
  }
  
  return surfaceArea / cbrt(36.0 * M_PI * totalVolume * totalVolume);
}

double DoubleCubicLatticeVolume::getPackingDensity() {

  if (vdwVolume == 0.0){
    getVDWVolume();
  }

  if (totalVolume == 0.0){
    getVolume();
  }

  return vdwVolume / totalVolume;
}

} // namespace Descriptors
} // namespace RDKit