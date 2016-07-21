//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../../RDGeneral/types.h"
#include "../../Geometry/point.h"
#include "Utilites.h"

namespace RDKit {
    namespace StructureCheck {

        void SetupNeighbourhood(const ROMol &mol, std::vector<Neighbourhood> &neighbours) {
            neighbours.clear();
            neighbours.resize(mol.getNumAtoms());

            for (unsigned i = 0; i < mol.getNumBonds(); i++) {
                const Bond* bond = mol.getBondWithIdx(i);
                unsigned a1 = bond->getBeginAtomIdx();
                unsigned a2 = bond->getEndAtomIdx();

                neighbours[a1].Atoms.push_back(a2);
                neighbours[a1].Bonds.push_back(i);

                neighbours[a2].Atoms.push_back(a1);
                neighbours[a2].Bonds.push_back(i);
            }
        }

        bool getMolAtomPoints(const ROMol& mol, std::vector<RDGeom::Point3D> &atomPoint) {
            bool non_zero_z = false;
            atomPoint.resize(mol.getNumAtoms());
            // take X,Y,Z coordinates of each atom
            if (0 != mol.getNumConformers())
                for (RDKit::ROMol::ConstConformerIterator cnfi = mol.beginConformers(); cnfi != mol.endConformers(); cnfi++) {
                    const Conformer &conf = **cnfi;  //mol.getConformer(confId);
                    if (conf.is3D()) {
                        for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
                            atomPoint[i] = conf.getAtomPos(i);
                            if (fabs(atomPoint[i].z) >= 1.e-7)
                                non_zero_z = true;
                        }
                        break;
                    }
                }
            if (atomPoint.empty()) {    // compute XYZ
//TODO:
// ???? ..........
            }
            return non_zero_z;
        }
 }// namespace StructureCheck
} // namespace RDKit
