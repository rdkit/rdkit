//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../MolOps.h"
#include "../Descriptors/MolDescriptors.h"
#include "StripSmallFragments.h"

namespace RDKit {
 namespace StructureCheck {

    static inline 
    std::string getMolecularFormula(const ROMol &mol) {
         return RDKit::Descriptors::calcMolFormula(mol);
    }

    void AddMWMF(RWMol &mol, bool pre) { // set formula & mass properties "MW_PRE" "MW_POST"
         double mass = 0.0;
         mass = RDKit::Descriptors::calcExactMW(mol);
/*
        for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
             const Atom& atom = *mol.getAtomWithIdx(i);
             mass += atom.getMass();
             mass += atom.getNumImplicitHs() * 1.0080; // and add implicit Hydrogens mass
         }
*/
         std::string formula = getMolecularFormula(mol);
         if(!formula.empty())
            mol.setProp((pre ? "MF_PRE" : "MF_POST"), formula);
         char propertyValue[64];
         snprintf(propertyValue, sizeof(propertyValue), "%g", mass);
         mol.setProp((pre ? "MW_PRE" : "MW_POST"), mass);
    }

     bool StripSmallFragments(RWMol &mol) {
        bool removed = false;

        std::vector<int>                        frags;
        std::vector<std::vector<int> >          fragsMolAtomMapping;
        std::vector<boost::shared_ptr<ROMol> >  molFrags = 
            RDKit::MolOps::getMolFrags(mol, true, &frags, &fragsMolAtomMapping);
        // remove all fragments smaller than the biggest one
        unsigned nf = fragsMolAtomMapping.size(); // for each fragment: RDKit::MolOps::getMolFrags(mol, frags);
        if (nf > 1) {
            unsigned maxFragSize = 0;
            unsigned maxSizeFragIdx = 0;
            for (unsigned i = 0; i < nf; i++) {
                if (molFrags[i]->getNumAtoms() > maxFragSize) {
                    maxFragSize = molFrags[i]->getNumAtoms();
                    maxSizeFragIdx = i;
                }
                for (unsigned i = 0; i < molFrags.size(); i++)
                    if(i != maxSizeFragIdx) {
                        for (unsigned j = 0; j < molFrags[i]->getNumAtoms(); i++)
                            mol.removeAtom(j);
                    }
            }
            removed = true;
        }
        return removed;
    }

 }// namespace StructureCheck
} // namespace RDKit
