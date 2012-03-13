// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "GasteigerCharges.h"
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include "GasteigerParams.h"

namespace Gasteiger {
  using namespace RDKit;
  /*! \brief split the formal charge across atoms of same type if we have a conjugated system
   *
   *  This function is called before the charge equivalization iteration is started for the 
   *  Gasteiger charges. If any of the atom involved in conjugated system have formal charges 
   *  set on them, this charges is equally distributed across atoms of the same type in that
   *  conjugated system. So for example the two nitrogens in the benzamidine system start the iteration
   * with equal charges of 0.5
   */
  void splitChargeConjugated(const ROMol &mol, DOUBLE_VECT &charges) {
    int aix;
    int natms = mol.getNumAtoms();
    INT_VECT marker;
    INT_VECT_CI mci;
    int aax, yax;
    double formal;
    const Atom *at, *aat, *yat;
    for (aix = 0; aix < natms; aix++) {
      at = mol.getAtomWithIdx(aix);
      formal = at->getFormalCharge();
      //std::cout << aix << " formal charges:" << formal << "\n";
      marker.resize(0);
      if ((fabs(formal) > EPS_DOUBLE)  && (fabs(charges[aix]) < EPS_DOUBLE) ) {
        marker.push_back(aix);
        ROMol::OEDGE_ITER bnd1, end1, bnd2, end2;
        boost::tie(bnd1,end1) = mol.getAtomBonds(at);
        while (bnd1 != end1) {
          if (mol[*bnd1]->getIsConjugated()) {
            aax = mol[*bnd1]->getOtherAtomIdx(aix);
            aat = mol.getAtomWithIdx(aax);
            boost::tie(bnd2,end2) = mol.getAtomBonds(aat);
            while (bnd2 != end2) {
              if ((*bnd1) != (*bnd2)) {
                if (mol[*bnd2]->getIsConjugated()) {
                  yax = mol[*bnd2]->getOtherAtomIdx(aax);
                  yat = mol.getAtomWithIdx(yax);
                  if (at->getAtomicNum() == yat->getAtomicNum()) {
                    formal += yat->getFormalCharge();
                    marker.push_back(yax);
                  }
                }
              }
              bnd2++;
            }
          }
          bnd1++;
        }

        for (mci = marker.begin(); mci != marker.end(); mci++) {
          charges[*mci] = (formal/marker.size());
        }
      }
    }
    /*
    for (aix = 0; aix < natms; aix++) {
      std::cout << "In splitter: " << " charges:" << charges[aix] << "\n";
      }*/
  }            
} // end of namespace Gasteiger

namespace RDKit {
  void computeGasteigerCharges(const ROMol *mol,int nIter, bool throwOnParamFailure) {
    PRECONDITION(mol,"bad molecule");
    computeGasteigerCharges(*mol,nIter,throwOnParamFailure);
  }
  void computeGasteigerCharges(const ROMol &mol,int nIter, bool throwOnParamFailure) {
    std::vector<double> chgs(mol.getNumAtoms());
    computeGasteigerCharges(mol,chgs,nIter,throwOnParamFailure);
  }


  /*! \brief compute the Gasteiger partial charges and return a new molecule with the charges set
   *
   * Ref : J.Gasteiger, M. Marsili, "Iterative Equalization of Oribital Electronegatiity
   *  A Rapid Access to Atomic Charges", Tetrahedron Vol 36 p3219 1980
   */
  void computeGasteigerCharges(const ROMol &mol, std::vector<double> &charges,
                               int nIter, bool throwOnParamFailure) {
    PRECONDITION(charges.size()>=mol.getNumAtoms(),"bad array size");
    
    PeriodicTable *table = PeriodicTable::getTable();
    const GasteigerParams *params = GasteigerParams::getParams();

    double damp = DAMP;
    int natms = mol.getNumAtoms();
    // space for parameters for each atom in the molecule 
    std::vector<DOUBLE_VECT> atmPs;
    atmPs.reserve(natms);
    
    std::fill(charges.begin(),charges.end(),0.0);

    DOUBLE_VECT hChrg; // total charge on the implicit hydrogen on each heavy atom
    hChrg.resize(natms, 0.0);

    DOUBLE_VECT ionX;
    ionX.resize(natms, 0.0);

    DOUBLE_VECT energ;
    energ.resize(natms, 0.0);

    ROMol::ADJ_ITER nbrIdx,endIdx;

    // deal with the conjugated system - distribute the formal charges on atoms of same type in each
    // conjugated system
    Gasteiger::splitChargeConjugated(mol, charges);

    // now read in the parameters
    ROMol::ConstAtomIterator ai; 
    
    for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
      std::string elem = table->getElementSymbol((*ai)->getAtomicNum());
      std::string mode;

      switch ((*ai)->getHybridization()) {
      case Atom::SP3:
        mode = "sp3";
        break;
      case Atom::SP2:
        mode = "sp2";
        break;
      case Atom::SP:
        mode = "sp";
        break;
      default:
        if ((*ai)->getAtomicNum() == 1) {
          // if it is hydrogen 
          mode = "*";
        } else if ((*ai)->getAtomicNum() == 16) {
          // we have a sulfur atom with no hydribidation information
          // check how many oxygens we have on the sulfer
          boost::tie(nbrIdx,endIdx) = mol.getAtomNeighbors(*ai);
          int no = 0;
          while (nbrIdx != endIdx) {
            if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() == 8){
              no++;
            }
            nbrIdx++;
          }
          if (no == 2) {
            mode = "so2";
          } else if (no == 1) {
            mode = "so";
          }
        }
      }

      // if we get a unknown mode or element type the 
      // following will will throw an exception
      atmPs.push_back(params->getParams(elem, mode,throwOnParamFailure));

      // set ionX paramters
      // if Hydrogen treat differently 
      int idx = (*ai)->getIdx();
      if ((*ai)->getAtomicNum() == 1) {
        ionX[idx] = IONXH;
      }
      else {
        ionX[idx] = atmPs[idx][0] + atmPs[idx][1] + atmPs[idx][2];
      }
    }

    // do the iteration here
    int itx, aix, sgn, niHs;
    double enr, dq, dx, qHs, dqH;
    // parameters for hydrogen atoms (for case where the hydrogen are not in the
    // graph (implicit hydrogens)
    DOUBLE_VECT hParams;
    hParams = params->getParams("H", "*", throwOnParamFailure);

    /*
    int itmp;
    for (itmp = 0; itmp < 5; itmp++) {
      std::cout << " aq:" << charges[itmp] << "\n";
      }*/

    for (itx = 0; itx < nIter; itx++) {
      for (aix = 0; aix < natms; aix++) {
        //enr = p0 + charge*(p1 + p2*charge)
        enr = atmPs[aix][0] + charges[aix]*(atmPs[aix][1] + atmPs[aix][2]*charges[aix]);
        energ[aix] = enr;
      }
      
      for (aix = 0; aix < natms; aix++) {
        dq = 0.0;
        boost::tie(nbrIdx,endIdx) = mol.getAtomNeighbors(mol.getAtomWithIdx(aix));
        while (nbrIdx != endIdx) {
          dx = energ[*nbrIdx] - energ[aix];
          if (dx < 0.0) {
            sgn = 0;
          } else {
            sgn = 1;
          }
          dq += dx/( (sgn*(ionX[aix] - ionX[*nbrIdx])) + ionX[*nbrIdx]);            
          nbrIdx++;
        }
        // now loop over the implicit hydrogens and get their contributions
        // since hydrogens don't connect to anything else, update their charges at the same time
        niHs = mol.getAtomWithIdx(aix)->getTotalNumHs();
        if (niHs > 0) {
          qHs = hChrg[aix]/niHs;
          enr = hParams[0] + qHs*(hParams[1] + hParams[2]*qHs);
          dx = enr - energ[aix];
          if (dx < 0.0) {
            sgn = 0;
          } else {
            sgn = 1;
          }
          
          dqH = dx/((sgn*(ionX[aix] - IONXH)) + IONXH);
          
          dq += (niHs*dqH);
          
          //adjust the charges on the hydrogens simultaneously (possible because each of the
          // hydrogens have no other neighbors)
          hChrg[aix] -= (niHs*dqH*damp);
        }
        charges[aix] += (damp*dq);
        
      }

      damp *= DAMP_SCALE;
    }

    for (aix = 0; aix < natms; aix++) {
      mol.getAtomWithIdx(aix)->setProp("_GasteigerCharge", charges[aix], true);
      // set the implicit hydrogen charges
      mol.getAtomWithIdx(aix)->setProp("_GasteigerHCharge", hChrg[aix], true);
      
    }
  }
}
          
  
