// $Id$
//
//  Copyright (C) 2013-2014 Paolo Tosco
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
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
  namespace MolAlign {
    typedef struct O3AFuncData {
      const Conformer *prbConf;
      const Conformer *refConf;
      void *prbProp;
      void *refProp;
      int coeff;
      int weight;
      bool useMMFFSim;
    } O3AFuncData;
    inline bool isDoubleZero(const double x) {
      return ((x < 1.0e-10) && (x > -1.0e-10));
    };
    
    class O3AConstraintVect;
    
    //! A class to define alignment constraints. Each constraint
    //! is defined by a pair of atom indexes (one for the probe,
    //! one for the reference) and a weight. Constraints can
    //! can be added via the O3AConstraintVect class.
    class O3AConstraint {
      friend class O3AConstraintVect;
      public:
        double getPrbIdx() {
          return d_prbIdx;
        }
        double getRefIdx() {
          return d_refIdx;
        }
        double getWeight() {
          return d_weight;
        }
      private:
        unsigned int d_idx;
        unsigned int d_prbIdx;
        unsigned int d_refIdx;
        double d_weight;
    };

    //! A class to store a vector of alignment constraints. Each constraint
    //! is defined by an O3AConstraint object. Each time the append()
    //! method is invoked, the vector is sorted to make lookup faster.
    //! Hence, constraints are not necessarily stored in the same order
    //! they were appended.
    class O3AConstraintVect {
      public:
        O3AConstraintVect() :
          d_count(0) {};
        ~O3AConstraintVect() {};
        void append(unsigned int prbIdx, unsigned int refIdx, double weight) {
          O3AConstraint *o3aConstraint = new O3AConstraint();
          o3aConstraint->d_idx = d_count;
          o3aConstraint->d_prbIdx = prbIdx;
          o3aConstraint->d_refIdx = refIdx;
          o3aConstraint->d_weight = weight;
          d_o3aConstraintVect.push_back(boost::shared_ptr<O3AConstraint>(o3aConstraint));
          std::sort(d_o3aConstraintVect.begin(),
            d_o3aConstraintVect.end(), d_compareO3AConstraint);
          ++d_count;
        }
        std::vector<boost::shared_ptr<O3AConstraint> >::size_type size() {
          return d_o3aConstraintVect.size();
        }
        O3AConstraint *operator[](unsigned int i) {
          return d_o3aConstraintVect[i].get();
        }
      private:
        unsigned int d_count;
        std::vector<boost::shared_ptr<O3AConstraint> > d_o3aConstraintVect;
        static bool d_compareO3AConstraint(boost::shared_ptr<O3AConstraint> a, boost::shared_ptr<O3AConstraint> b)
        {
          return ((a->d_prbIdx != b->d_prbIdx) ? (a->d_prbIdx < b->d_prbIdx)
            : ((a->d_refIdx != b->d_refIdx) ? (a->d_refIdx < b->d_refIdx)
            : ((a->d_weight != b->d_weight) ? (a->d_weight > b->d_weight)
            : (a->d_idx < b->d_idx))));
        };
    };

    const int O3_DUMMY_COST = 100000;
    const int O3_LARGE_NEGATIVE_WEIGHT = -1e9;
    const int O3_DEFAULT_CONSTRAINT_WEIGHT = 100.0;
    const unsigned int O3_MAX_H_BINS = 20;
    const unsigned int O3_MAX_SDM_ITERATIONS = 100;
    const unsigned int O3_MAX_SDM_THRESHOLD_ITER = 3;
    const double O3_RANDOM_TRANS_COEFF = 5.0;
    const double O3_THRESHOLD_DIFF_DISTANCE = 0.1;
    const double O3_SDM_THRESHOLD_START = 0.7;
    const double O3_SDM_THRESHOLD_STEP = 0.3;
    const double O3_CHARGE_WEIGHT = 10.0;
    const double O3_CRIPPEN_WEIGHT = 10.0;
    const double O3_RMSD_THRESHOLD = 1.0e-04;
    const double O3_SCORE_THRESHOLD = 0.01;
    const double O3_SCORING_FUNCTION_ALPHA = 5.0;
    const double O3_SCORING_FUNCTION_BETA = 0.5;
    const double O3_CHARGE_COEFF = 5.0;
    const double O3_CRIPPEN_COEFF = 1.0;
    const int O3_MAX_WEIGHT_COEFF = 5;
    enum {
      O3_USE_MMFF_WEIGHTS = (1<<0),
      O3_ACCURACY_MASK = (1<<0 | 1<<1),
      O3_LOCAL_ONLY = (1<<2)
    };
    
    class MolHistogram {
    public:
      MolHistogram(const ROMol &mol, const double *dmat);
      ~MolHistogram() {};
    inline int get(const unsigned int y, const unsigned int x) const {
      PRECONDITION(y < d_h.shape()[0], "Invalid index on MolHistogram");
      PRECONDITION(x < d_h.shape()[1], "Invalid index on MolHistogram");
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
        PRECONDITION(i < d_cost.shape()[0], "Invalid index on LAP.cost");
        PRECONDITION(j < d_cost.shape()[1], "Invalid index on LAP.cost");
        return d_cost[i][j];
      }
      int getRowSol(const unsigned int i) {
        PRECONDITION(i < d_rowSol.size(), "Invalid index on LAP.rowSol");
        return d_rowSol[i];
      }
      void computeMinCostPath(const int dim);
      void computeCostMatrix(const ROMol &prbMol, const MolHistogram &prbHist,
        const ROMol &refMol, const MolHistogram &refHist,
        O3AConstraintVect *o3aConstraintVect, int (*costFunc)
        (const unsigned int, const unsigned int, double, void *),
        void *data, const unsigned int n_bins = O3_MAX_H_BINS);
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
      SDM(const Conformer *prbConf = NULL, const Conformer *refConf = NULL,
        O3AConstraintVect *o3aConstraintVect = NULL) : 
        d_prbConf(prbConf),
        d_refConf(refConf),
        d_o3aConstraintVect(o3aConstraintVect) {};
      // copy constructor
      SDM(const SDM &other) :
        d_prbConf(other.d_prbConf),
        d_refConf(other.d_refConf),
        d_o3aConstraintVect(other.d_o3aConstraintVect),
        d_SDMPtrVect(other.d_SDMPtrVect.size()) {
        for (unsigned int i = 0; i < d_SDMPtrVect.size(); ++i) {
          d_SDMPtrVect[i] = boost::shared_ptr<SDMElement>(new SDMElement());
          memcpy(d_SDMPtrVect[i].get(), other.d_SDMPtrVect[i].get(), sizeof(SDMElement));
        }
      };
      // assignment operator
      SDM& operator=(const SDM &other) {
        d_prbConf = other.d_prbConf;
        d_refConf = other.d_refConf;
        d_o3aConstraintVect = other.d_o3aConstraintVect;
        d_SDMPtrVect.resize(other.d_SDMPtrVect.size());
        for (unsigned int i = 0; i < d_SDMPtrVect.size(); ++i) {
          d_SDMPtrVect[i] = boost::shared_ptr<SDMElement>(new SDMElement());
          memcpy(d_SDMPtrVect[i].get(), other.d_SDMPtrVect[i].get(), sizeof(SDMElement));
        }
        
        return *this;
      };
      // destructor    
      ~SDM() {};
      void fillFromDist(double threshold,
                        const boost::dynamic_bitset<> &refHvyAtoms,
                        const boost::dynamic_bitset<> &prbHvyAtoms);
      void fillFromLAP(LAP &lap);
      double scoreAlignment(double (*scoringFunc)
        (const unsigned int, const unsigned int, void *), void *data);
      void prepareMatchWeightsVect(RDKit::MatchVectType &matchVect,
        RDNumeric::DoubleVector &weights, double (*weightFunc)
        (const unsigned int, const unsigned int, void *), void *data);
      unsigned int size() {
        return d_SDMPtrVect.size();
      }
   private:
      typedef struct SDMElement {
        unsigned int idx[2];
        int score;
        int cost;
        double sqDist;
        O3AConstraint *o3aConstraint;
      } SDMElement;
      const Conformer *d_prbConf;
      const Conformer *d_refConf;
      O3AConstraintVect *d_o3aConstraintVect;
      std::vector<boost::shared_ptr<SDMElement> > d_SDMPtrVect;
      static bool compareSDMScore(boost::shared_ptr<SDMElement> a, boost::shared_ptr<SDMElement> b)
      {
        return ((a->score != b->score) ? (a->score < b->score)
          : ((a->cost != b->cost) ? (a->cost < b->cost)
          : ((a->idx[0] != b->idx[0]) ? (a->idx[0] < b->idx[0])
          : (a->idx[1] < b->idx[1]))));
      };
      static bool compareSDMDist(boost::shared_ptr<SDMElement> a, boost::shared_ptr<SDMElement> b)
      {
        double aWeight = (a->o3aConstraint ? a->o3aConstraint->getWeight() : 0.0);
        double bWeight = (b->o3aConstraint ? b->o3aConstraint->getWeight() : 0.0);
        return ((aWeight != bWeight) ? (aWeight > bWeight)
          : ((a->sqDist != b->sqDist) ? (a->sqDist < b->sqDist)
          : ((a->idx[0] != b->idx[0]) ? (a->idx[0] < b->idx[0])
          : (a->idx[1] < b->idx[1]))));
      };
    };
    
    class O3A {
    public:
      //! pre-defined atom typing schemes
      typedef enum {
        MMFF94=0,
        CRIPPEN
      } AtomTypeScheme;
      O3A(ROMol &prbMol, const ROMol &refMol,
          void *prbProp, void *refProp, AtomTypeScheme atomTypes = MMFF94,
          const int prbCid = -1, const int refCid = -1,
          const bool reflect = false, const unsigned int maxIters = 50,
          unsigned int options = 0, const MatchVectType *constraintMap = NULL,
          const RDNumeric::DoubleVector *constraintWeights = NULL, LAP *extLAP = NULL,
          MolHistogram *extPrbHist = NULL, MolHistogram *extRefHist = NULL);
      O3A(int (*costFunc)(const unsigned int, const unsigned int, double, void *),
        double (*weightFunc)(const unsigned int, const unsigned int, void *),
        double (*scoringFunc)(const unsigned int, const unsigned int, void *),
        void *data, ROMol &prbMol, const ROMol &refMol,
        const int prbCid, const int refCid,
        boost::dynamic_bitset<> *prbHvyAtoms = NULL,
        boost::dynamic_bitset<> *refHvyAtoms = NULL,
        const bool reflect = false, const unsigned int maxIters = 50,
        unsigned int options = 0, O3AConstraintVect *o3aConstraintVect = NULL,
        ROMol *extWorkPrbMol = NULL, LAP *extLAP = NULL,
        MolHistogram *extPrbHist = NULL, MolHistogram *extRefHist = NULL);
      ~O3A() {
        if (d_o3aMatchVect) {
          delete d_o3aMatchVect;
        }
        if (d_o3aWeights) {
          delete d_o3aWeights;
        }
      };
      double align();
      double trans(RDGeom::Transform3D &trans);
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
    int o3aMMFFCostFunc(const unsigned int prbIdx, const unsigned int refIdx, double hSum, void *data);
    double o3aMMFFWeightFunc(const unsigned int prbIdx, const unsigned int refIdx, void *data);
    double o3aMMFFScoringFunc(const unsigned int prbIdx, const unsigned int refIdx, void *data);
    int o3aCrippenCostFunc(const unsigned int prbIdx, const unsigned int refIdx, double hSum, void *data);
    double o3aCrippenWeightFunc(const unsigned int prbIdx, const unsigned int refIdx, void *data);
    double o3aCrippenScoringFunc(const unsigned int prbIdx, const unsigned int refIdx, void *data);
  }
}
#endif
