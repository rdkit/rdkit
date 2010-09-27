// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SHAPE_UTILS_H_20050128_
#define _RD_SHAPE_UTILS_H_20050128_
#include <DataStructs/DiscreteValueVect.h>
#include <vector>

namespace RDGeom {
  class Point3D;
  class Transform3D;
}

namespace RDKit {
  class ROMol;
  class Conformer;

  namespace MolShapes {

    //! Compute the size of the box that can fit the conformation, and offset of the box
    //! from the origin
    void computeConfDimsAndOffset(const Conformer &conf, RDGeom::Point3D &dims, 
                              RDGeom::Point3D &offSet, const RDGeom::Transform3D *trans=0,
                              double padding=2.5);
    
    //! Compute the a box that will fit the confomer
    /*!
      \param conf            The conformer of interest
      \param leftBottom      Storage for one extremity of the box
      \param rightTop        Storage for other extremity of the box
      \param trans           Optional transformation to be applied to the atom coordinates
      \param padding         Padding added on the sides around the conformer
    */
    void computeConfBox(const Conformer &conf, RDGeom::Point3D &leftBottom, 
                         RDGeom::Point3D &rightTop, const RDGeom::Transform3D *trans=0,
                        double padding=2.5);

    //! Compute the union of two boxes
    void computeUnionBox(const RDGeom::Point3D &leftBottom1, const RDGeom::Point3D &rightTop1, 
                         const RDGeom::Point3D &leftBottom2, const RDGeom::Point3D &rightTop2, 
                         RDGeom::Point3D &uLeftBottom, RDGeom::Point3D &uRightTop);

    //! Compute dimensions of a conformer
    /*!
      \param conf     Conformer of interest
      \param padding  Padding added to the atom coordinates on all sides
      \param center   Optionally specify the center 
      \param ignoreHs if true, ignore the hydrogen atoms in computing the centroid
    */
    std::vector<double> getConfDimensions(const Conformer &conf, double padding=2.5, 
                                          const RDGeom::Point3D *center=0, bool ignoreHs=true);

    //! Compute the shape tanimoto distance between two molecule based on a predefined alignment
    /*!
      \param mol1         The first molecule of interest
      \param mol2         The second molecule of interest
      \param confId1      Conformer in the first molecule (defaults to first conformer)
      \param confId2      Conformer in the second molecule (defaults to first conformer)
      \param gridSpacing  resolution of the grid used to encode the molecular shapes
      \param bitsPerPoint number of bit used to encode the occupancy at each grid point
                          defaults to two bits per grid point
      \param vdwScale     Scaling factor for the radius of the atoms to determine the base radius 
                          used in the encoding - grid points inside this sphere carry the maximum occupany
      \param stepSize     thickness of the each layer outside the base radius, the occupancy value is decreased 
                          from layer to layer from the maximum value
      \param maxLayers    the maximum number of layers - defaults to the number allowed the number of bits 
                          use per grid point - e.g. two bits per grid point will allow 3 layers
      \param ignoreHs     if true, ignore the hydrogen atoms in the shape encoding process
     */

    double tanimotoDistance(const ROMol &mol1, const ROMol &mol2, int confId1=-1, int confId2=-1,
                            double gridSpacing=0.5, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE, 
                            double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1,
                            bool ignoreHs=true);


    //! Compute the shape tanimoto distance between two conformers based on a predefined alignment
    /*!
      \param conf1        The first conformer of interest
      \param conf2        The second conformer of interest
      \param gridSpacing  resolution of the grid used to encode the molecular shapes
      \param bitsPerPoint number of bit used to encode the occupancy at each grid point
      \param vdwScale     Scaling factor for the radius of the atoms to determine the base radius 
                          used in the encoding - grid points inside this sphere carry the maximum occupany
      \param stepSize     thickness of the each layer outside the base radius, the occupancy value is decreased 
                          from layer to layer from the maximum value
      \param maxLayers    the maximum number of layers - defaults to the number allowed the number of bits 
                          use per grid point - e.g. two bits per grid point will allow 3 layers
      \param ignoreHs     if true, ignore the hydrogen atoms in the shape encoding process
     */

    double tanimotoDistance(const Conformer &conf1, const Conformer &conf2, double gridSpacing=0.5, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE,
                            double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1, bool ignoreHs=true);


    //! Compute the shape protrusion distance between two molecule based on a predefined alignment
    /*!
      \param mol1         The first molecule of interest
      \param mol2         The second molecule of interest
      \param confId1      Conformer in the first molecule (defaults to first conformer)
      \param confId2      Conformer in the second molecule (defaults to first conformer)
      \param gridSpacing  resolution of the grid used to encode the molecular shapes
      \param bitsPerPoint number of bit used to encode the occupancy at each grid point
                          defaults to two bits per grid point
      \param vdwScale     Scaling factor for the radius of the atoms to determine the base radius 
                          used in the encoding - grid points inside this sphere carry the maximum occupany
      \param stepSize     thickness of the each layer outside the base radius, the occupancy value is decreased 
                          from layer to layer from the maximum value
      \param maxLayers    the maximum number of layers - defaults to the number allowed the number of bits 
                          use per grid point - e.g. two bits per grid point will allow 3 layers
      \param ignoreHs     if true, ignore the hydrogen atoms in the shape encoding process
      \param allowReordering  if set the order will be automatically updated so that the value calculated
                              is the protrusion of the smaller shape from the larger one.
     */

    double protrudeDistance(const ROMol &mol1, const ROMol &mol2, int confId1=-1, int confId2=-1,
                            double gridSpacing=0.5, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE, 
                            double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1,
                            bool ignoreHs=true, bool allowReordering=true);


    //! Compute the shape protrusion distance between two conformers based on a predefined alignment
    /*!
      \param conf1        The first conformer of interest
      \param conf2        The second conformer of interest
      \param gridSpacing  resolution of the grid used to encode the molecular shapes
      \param bitsPerPoint number of bit used to encode the occupancy at each grid point
      \param vdwScale     Scaling factor for the radius of the atoms to determine the base radius 
                          used in the encoding - grid points inside this sphere carry the maximum occupany
      \param stepSize     thickness of the each layer outside the base radius, the occupancy value is decreased 
                          from layer to layer from the maximum value
      \param maxLayers    the maximum number of layers - defaults to the number allowed the number of bits 
                          use per grid point - e.g. two bits per grid point will allow 3 layers
      \param ignoreHs     if true, ignore the hydrogen atoms in the shape encoding process
      \param allowReordering  if set the order will be automatically updated so that the value calculated
                              is the protrusion of the smaller shape from the larger one.
     */

    double protrudeDistance(const Conformer &conf1, const Conformer &conf2, double gridSpacing=0.5, 
                            DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE,
                            double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1, bool ignoreHs=true,
                            bool allowReordering=true);
    
  }

}

#endif
