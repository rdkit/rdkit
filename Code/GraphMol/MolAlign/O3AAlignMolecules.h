//
//  Copyright (C) 2013 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_O3AALIGNMOLECULES_H_
#define _RD_O3AALIGNMOLECULES_H_

#include <RDGeneral/Invariant.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <Numerics/Vector.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <vector>
#include <cmath>
#include "boost/multi_array.hpp"

namespace RDKit {
  namespace MolAlign {
    inline const bool isDoubleZero(const double x) {
      return ((x < 1.0e-10) && (x > -1.0e-10));
    };
    
    const int O3_DUMMY_COST = 100000;
    const unsigned int O3_MAX_H_BINS = 20;
    const unsigned int O3_MAX_SDM_ITERATIONS = 100;
    const double O3_RANDOM_TRANS_COEFF = 5.0;
    const double O3_THRESHOLD_DIFF_DISTANCE = 0.1;
    const double O3_SDM_THRESHOLD = 0.7;
    const double O3_CHARGE_WEIGHT = 10.0;
    const double O3_RMSD_THRESHOLD = 1.0e-04;
    const double O3_SCORE_THRESHOLD = 1.0e-04;
    const double O3_SCORING_FUNCTION_ALPHA = 5.0;
    const double O3_SCORING_FUNCTION_BETA = 0.5;
    const double O3_CHARGE_COEFF = 5.0;
    const int O3_MAX_WEIGHT_COEFF = 5;
    enum {
      O3_USE_MMFF_WEIGHTS = (1<<0)
    };
    
    class MolHistogram {
    public:
      MolHistogram(const ROMol &mol, const int cid = -1);
      ~MolHistogram() {};
    inline int get(const unsigned int y, const unsigned int x) const {
      PRECONDITION((y >= 0) && (y < d_h.shape()[0]), "Invalid index on MolHistogram");
      PRECONDITION((x >= 0) && (x < d_h.shape()[1]), "Invalid index on MolHistogram");
      return d_h[y][x];
    }
    private:
      boost::multi_array<int, 2> d_h;
    };
    
    class LAP {
    public:
      LAP(unsigned int dim) :
        d_rowSol(dim),
        d_colSol(dim),
        d_free(dim),
        d_colList(dim),
        d_matches(dim),
        d_d(dim),
        d_v(dim),
        d_pred(dim),
        d_cost(boost::extents[dim][dim]) {};
      ~LAP() {};
      int getCost(const unsigned int i, const unsigned int j) {
        PRECONDITION((i >= 0) && (i < d_cost.shape()[0]), "Invalid index on LAP.cost");
        PRECONDITION((j >= 0) && (j < d_cost.shape()[1]), "Invalid index on LAP.cost");
        return d_cost[i][j];
      }
      int getRowSol(const unsigned int i) {
        PRECONDITION((i >= 0) && (i < d_rowSol.size()), "Invalid index on LAP.rowSol");
        return d_rowSol[i];
      }
      void computeMinCostPath(unsigned int dim);
      void computeCostMatrix(const ROMol &prbMol, const MolHistogram &prbHist,
        MMFF::MMFFMolProperties *prbMP, const ROMol &refMol, const MolHistogram &refHist,
        MMFF::MMFFMolProperties *refMP, const int coeff = 0, const int n_bins = O3_MAX_H_BINS);
    private:
      std::vector<int> d_rowSol;
      std::vector<int> d_colSol;
      std::vector<int> d_free;
      std::vector<int> d_colList;
      std::vector<int> d_matches;
      std::vector<int> d_d;
      std::vector<int> d_v;
      std::vector<int> d_pred;
      boost::multi_array<int, 2> d_cost;
    };
    
    class SDM {
    public:
      // constructor
      SDM(const ROMol *prbMol = NULL, const ROMol *refMol = NULL,
        MMFF::MMFFMolProperties *prbMP = NULL, MMFF::MMFFMolProperties *refMP = NULL,
        const int prbCid = -1, const int refCid = -1, const bool reflect = false) : 
        d_prbMol(prbMol),
        d_refMol(refMol),
        d_prbMMFFMolProperties(prbMP),
        d_refMMFFMolProperties(refMP),
        d_prbCid(prbCid),
        d_refCid(refCid) {};
      // copy constructor
      SDM(const SDM &other) :
        d_prbMol(other.d_prbMol),
        d_refMol(other.d_refMol),
        d_prbMMFFMolProperties(other.d_prbMMFFMolProperties),
        d_refMMFFMolProperties(other.d_refMMFFMolProperties),
        d_prbCid(other.d_prbCid),
        d_refCid(other.d_refCid),
        d_SDMPtrVect(other.d_SDMPtrVect.size()) {
        for (unsigned int i = 0; i < d_SDMPtrVect.size(); ++i) {
          d_SDMPtrVect[i] = new SDMElement();
          memcpy(d_SDMPtrVect[i], other.d_SDMPtrVect[i], sizeof(SDMElement));
        }
      };
      // assignment operator
      SDM& operator=(const SDM &other) {
        d_prbMol = other.d_prbMol;
        d_refMol = other.d_refMol;
        d_prbMMFFMolProperties = other.d_prbMMFFMolProperties;
        d_refMMFFMolProperties = other.d_refMMFFMolProperties;
        d_prbCid = other.d_prbCid;
        d_refCid = other.d_refCid;
        d_SDMPtrVect.resize(other.d_SDMPtrVect.size());
        for (unsigned int i = 0; i < d_SDMPtrVect.size(); ++i) {
          d_SDMPtrVect[i] = new SDMElement();
          memcpy(d_SDMPtrVect[i], other.d_SDMPtrVect[i], sizeof(SDMElement));
        }
        
        return *this;
      };
      // destructor    
      ~SDM() {
        for (unsigned int i = 0; i < d_SDMPtrVect.size(); ++i) {
          delete d_SDMPtrVect[i];
        }
      };
      void fillFromDist(const double threshold = O3_SDM_THRESHOLD);
      void fillFromLAP(LAP &lap);
      double scoreAlignment();
      void prepareMatchWeightsVect(RDKit::MatchVectType &matchVect,
        RDNumeric::DoubleVector &weights, int coeff = 0);
      unsigned int size() {
        return d_SDMPtrVect.size();
      }
   private:
      typedef struct SDMElement {
        unsigned int idx[2];
        int score;
        int cost;
        double sqDist;
      } SDMElement;
      const ROMol *d_prbMol;
      const ROMol *d_refMol;
      MMFF::MMFFMolProperties *d_prbMMFFMolProperties;
      MMFF::MMFFMolProperties *d_refMMFFMolProperties;
      int d_prbCid;
      int d_refCid;
      std::vector<SDMElement *> d_SDMPtrVect;
      static bool compareSDMScore(SDMElement *a, SDMElement *b)
      {
        return ((a->score != b->score) ? (a->score < b->score)
          : ((a->cost != b->cost) ? (a->cost < b->cost)
          : (a->idx[0] < b->idx[0])));
      };
      static bool compareSDMDist(SDMElement *a, SDMElement *b)
      {
        return (isDoubleZero(a->sqDist - b->sqDist)
          ? (a->idx[0] < b->idx[0]) : (a->sqDist < b->sqDist));
      };
    };
    
    class O3A {
    public:
      O3A(ROMol &prbMol, const ROMol &refMol,
        MMFF::MMFFMolProperties *prbMP, MMFF::MMFFMolProperties *refMP,
        const int prbCid = -1, const int refCid = -1,
        const bool reflect = false, const unsigned int maxIters = 50,
        const unsigned int options = O3_USE_MMFF_WEIGHTS,
        LAP *extLAP = NULL);
      ~O3A() {
        if (d_o3aMatchVect) {
          delete d_o3aMatchVect;
        }
        if (d_o3aWeights) {
          delete d_o3aWeights;
        }
      };
      double align() {
        return alignMol(*d_prbMol, *d_refMol, d_prbCid, d_refCid,
          d_o3aMatchVect, d_o3aWeights, d_reflect, d_maxIters);
      };
      std::pair<double, RDGeom::Transform3D *> trans() {
        RDGeom::Transform3D *trans = new RDGeom::Transform3D();
        double rmsd = getAlignmentTransform(*d_prbMol, *d_refMol,
          *trans, d_prbCid, d_refCid, d_o3aMatchVect, d_o3aWeights,
          d_reflect, d_maxIters);
        
        return std::make_pair(rmsd, trans);
      };
      double score() {
        return d_o3aScore;
      };
      const RDKit::MatchVectType *matches() {
        return d_o3aMatchVect;
      };
      const RDNumeric::DoubleVector *weights() {
        return d_o3aWeights;
      };
    private:
      ROMol *d_prbMol;
      const ROMol *d_refMol;
      MMFF::MMFFMolProperties *d_prbMMFFMolProperties;
      MMFF::MMFFMolProperties *d_refMMFFMolProperties;
      int d_prbCid;
      int d_refCid;
      bool d_reflect;
      unsigned int d_maxIters;
      const RDKit::MatchVectType *d_o3aMatchVect;
      const RDNumeric::DoubleVector *d_o3aWeights;
      double d_o3aScore;
    };
    
    void randomTransform(ROMol &mol, const int cid = -1, const int seed = -1);
    const RDGeom::POINT3D_VECT *reflect(const Conformer &conf);
  }
}
#endif
