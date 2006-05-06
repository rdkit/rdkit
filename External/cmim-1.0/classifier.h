
////////////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or              //
// modify it under the terms of the GNU General Public License                //
// version 2 as published by the Free Software Foundation.                    //
//                                                                            //
// This program is distributed in the hope that it will be useful, but        //
// WITHOUT ANY WARRANTY; without even the implied warranty of                 //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          //
// General Public License for more details.                                   //
//                                                                            //
// Written by François Fleuret                                                //
// Contact <francois.fleuret@epfl.ch> for comments & bug reports              //
// Copyright (C) 2004 EPFL                                                    //
////////////////////////////////////////////////////////////////////////////////

// $Id: classifier.h,v 1.1 2005/03/03 15:52:35 fleuret Exp $

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

using namespace std;

#include <iostream>
#include <fstream>

#include "misc.h"
#include "fastentropy.h"

class FeatureSelector;

class DataSet {
  // This keeps track of the number of references
  struct RawData {
    int nrefs;
    uint32_t *x, *y;
  };
public:
  int nb_samples, nb_features;
  int size;
  RawData *raw;
  uint32_t *y_va, **x_va;
  DataSet(istream &is,bool delimited=true);
  DataSet(int nb_samples, int nb_features);
  DataSet(const DataSet &ds);
  DataSet(const DataSet &ds, const FeatureSelector &fs);
  DataSet(const DataSet &ds, bool *selected_samples);
  DataSet &operator = (const DataSet &ds);
  ~DataSet();
  void copy(const DataSet &ds);
  void save_ascii(ostream &os);
};

//////////////////////////////////////////////////////////////////////
// The classifier ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

class Classifier {
public:
  enum { ID_LINEAR };
  static Classifier *load(istream &is);
  virtual void predict(const DataSet &ds, float *result) = 0;
  virtual void save(ostream &out) = 0;
  virtual void inner_load(istream &in) = 0;
};

class FeatureSelector {
public:
  int nb_selected_features;
  int *selected_index;

  // Those remains from the feature selection process. They can be
  // used as-is in the case of adaboost
  float *weights;

  FeatureSelector(istream &is);
  FeatureSelector(int nb_selected_features);
  ~FeatureSelector();

  // Selects features according to the Conditional Mutual Information Maximisation
  void cmim(const DataSet &ds);

  // Selects features according to the Mutual Information Maximisation
  void mim(const DataSet &ds);

  // Selects random features
  void random(const DataSet &ds);

  void save(ostream &os);
};

class LinearClassifier : public Classifier {
  int nb_features;
  float *weights;
  float bias;
public:
  LinearClassifier();
  LinearClassifier(int nb_features);
  virtual ~LinearClassifier();

  void compute_bayesian_weights(int nb_samples, uint32_t *y_va, uint32_t **x_va);
  void compute_bias(int nb_samples, uint32_t *y_va, uint32_t **x_va, bool balanced_error);

  void learn_bayesian(const DataSet &ds, bool balanced_error);
  void learn_perceptron(const DataSet &ds, bool balanced_error);

  virtual void predict(const DataSet &ds, float *result);
  virtual void save(ostream &out);
  virtual void inner_load(istream &is);
};

void compute_error_rates(FeatureSelector *selector, Classifier *classifier,
                         const DataSet &testing_set, int &n00, int &n01, int &n10, int &n11, float *result);

#endif
