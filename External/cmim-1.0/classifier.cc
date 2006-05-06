
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

// $Id: classifier.cc,v 1.1 2005/03/03 15:52:35 fleuret Exp $

#include <stdlib.h>
#include <cmath>

#include "classifier.h"

DataSet::DataSet(const DataSet &ds) {
  copy(ds);
}

DataSet::DataSet(int ns, int nf) {
  nb_samples = ns;
  nb_features = nf;

  size = (nb_samples + 31)/32;

  raw = new RawData;

  raw->nrefs = 1;
  raw->x = new uint32_t[nb_features * size];
  raw->y = new uint32_t[size];

  y_va = raw->y;
  x_va = new uint32_t *[nb_features];
  for(int j = 0; j < nb_features; j++) x_va[j] = raw->x + size*j;
}

DataSet::DataSet(istream &is,bool delimited) {
  is >> nb_samples >> nb_features;
  //cerr << ">>>" << nb_samples << " " << nb_features << "\n"; cerr.flush();
  size = (nb_samples + 31)/32;

  raw = new RawData;

  raw->nrefs = 1;
  raw->x = new uint32_t[nb_features * size];
  raw->y = new uint32_t[size];

  y_va = raw->y;
  x_va = new uint32_t *[nb_features];
  for(int j = 0; j < nb_features; j++) x_va[j] = raw->x + size*j;

  int v;

  for(int s = 0; s < nb_samples; s++) {
    //if(!(s%500)) cerr << "s=" << s << "\n"; cerr.flush();
    for(int f = 0; f < nb_features; f++) {
      if(delimited){
	is >> v;
	fe_set_bit(s, x_va[f], v > 0);
      } else {
	char cV='\n';
	while(cV=='\n')	is.get(cV);
	fe_set_bit(s, x_va[f], cV=='1');
      }
    }

    is >> v;
    fe_set_bit(s, y_va, v > 0);

    if(is.eof()) {
      cerr << "Error: missing data!\n";
      exit(1);
    }
  }
}

DataSet::DataSet(const DataSet &ds, const FeatureSelector &fs) {
  nb_samples = ds.nb_samples;
  nb_features = fs.nb_selected_features;
  size = (nb_samples + 31)/32;

  raw = ds.raw;

  raw->nrefs++;

  y_va = raw->y;
  x_va = new uint32_t *[nb_features];
  for(int j = 0; j < nb_features; j++) x_va[j] = raw->x + size*fs.selected_index[j];
}

DataSet::DataSet(const DataSet &ds, bool *selected_samples) {
  nb_samples = 0;
  for(int i = 0; i < ds.nb_samples; i++) if(selected_samples[i]) nb_samples++;
  nb_features = ds.nb_features;
  size = (nb_samples + 31)/32;
  raw = new RawData;

  raw->nrefs = 1;
  raw->x = new uint32_t[nb_features * size];
  raw->y = new uint32_t[size];

  y_va = raw->y;
  x_va = new uint32_t *[nb_features];
  for(int j = 0; j < nb_features; j++) x_va[j] = raw->x + size*j;

  int k = 0;
  for(int s = 0; s < ds.nb_samples; s++) if(selected_samples[s]) {
    fe_set_bit(k, y_va, fe_get_bit(s, ds.y_va));
    for(int j = 0; j < nb_features; j++) fe_set_bit(k, x_va[j], fe_get_bit(s, ds.x_va[j]));
    k++;
  }
}

DataSet &DataSet::operator = (const DataSet &ds) {
  if(this != &ds) {
    delete[] x_va;
    raw->nrefs--;
    if(raw->nrefs == 0) {
      delete[] raw->x;
      delete[] raw->y;
      delete raw;
    }
    copy(ds);
  }
  return *this;
}

DataSet::~DataSet() {
  delete[] x_va;
  raw->nrefs--;
  if(raw->nrefs == 0) {
    delete[] raw->x;
    delete[] raw->y;
    delete raw;
  }
}

void DataSet::copy(const DataSet &ds) {
  nb_samples = ds.nb_samples;
  nb_features = ds.nb_features;
  size = (nb_samples + 31)/32;
  raw = ds.raw;
  raw->nrefs++;
  y_va = raw->y;
  x_va = new uint32_t *[nb_features];
  for(int j = 0; j < nb_features; j++) x_va[j] = ds.x_va[j];
}

void DataSet::save_ascii(ostream &os) {
  os << nb_samples << " " << nb_features << "\n";
  for(int s = 0; s < nb_samples; s++) {
    for(int f = 0; f < nb_features; f++) {
      if(fe_get_bit(s, x_va[f])) os << "1"; else os << "0";
      if(f < nb_features-1) os << " "; else os << "\n";
    }
    if(fe_get_bit(s, y_va)) os << "1\n"; else os << "0\n";
  }
}

Classifier *Classifier::load(istream &is) {
  Classifier *result;

  int scheme;
  is >> scheme;

  switch(scheme) {
  case Classifier::ID_LINEAR:
    result = new LinearClassifier();
    break;

  default:
    result = 0;
    cerr << "Unknown classifier type!\n";
    exit(1);
  }

  result->inner_load(is);

  return result;
}

FeatureSelector::FeatureSelector(int nb) : nb_selected_features(nb),
                                           selected_index(new int[nb_selected_features]),
                                           weights(new float[nb_selected_features]) { }

FeatureSelector::~FeatureSelector() {
  delete[] weights;
  delete[] selected_index;
}

LinearClassifier::LinearClassifier(int nf) : nb_features(nf),
                                             weights(new float[nb_features]),
                                             bias(0) { }

LinearClassifier::LinearClassifier() : nb_features(0),
                                       weights(0),
                                       bias(0) { }

LinearClassifier::~LinearClassifier() {
  delete[] weights;
}

void LinearClassifier::compute_bayesian_weights(int nb_samples, uint32_t *y_va, uint32_t **x_va) {
  for(int nf = 0; nf < nb_features; nf++) {
    int n11 = fe_count_and(nb_samples, y_va, x_va[nf]);
    int n10 = fe_count_and_not(nb_samples, y_va, x_va[nf]);
    int n01 = fe_count_and_not(nb_samples, x_va[nf], y_va);
    int n00 = fe_count_and_not_not(nb_samples, y_va, x_va[nf]);
    if(n00 == 0) n00 = 1; // This is sort of a dirty way to emulate +/- infty
    if(n01 == 0) n01 = 1;
    if(n10 == 0) n10 = 1;
    if(n11 == 0) n11 = 1;
    weights[nf] = log(float(n11 * n00) / float(n10 * n01));
  }
}

void LinearClassifier::compute_bias(int nb_samples, uint32_t *y_va, uint32_t **x_va, bool balanced_error) {
  Couple tmp[nb_samples];

  int n00 = 0, n01 = 0, n10 = 0, n11 = 0;
  for(int s = 0; s < nb_samples; s++) {
    float r = 0;
    for(int nf = 0; nf < nb_features; nf++) if(fe_get_bit(s, x_va[nf])) r += weights[nf];
    tmp[s].index = s;
    tmp[s].value = r;
    if(fe_get_bit(s, y_va)) n11++; else n01++;
  }

  qsort(tmp, nb_samples, sizeof(Couple), compare_couple);

  float error, best_error = 2.0;
  bias = 0;
  for(int t = 0; t < nb_samples-1; t++) {
    if(fe_get_bit(tmp[t].index, y_va)) { n10++; n11--; } else { n01--; n00++; }

    if(balanced_error)
      error = (float(n01) / float(n00 + n01) + float(n10) / float(n10 + n11)) / 2.0;
    else
      error = float(n01+n10) / float(n00 + n01 + n10 + n11);

    if(error < best_error) {
      best_error = error;
      bias = - (tmp[t].value + tmp[t+1].value)/2.0;
    }
  }
}

void LinearClassifier::learn_bayesian(const DataSet &ds, bool balanced_error) {
  compute_bayesian_weights(ds.nb_samples, ds.y_va, ds.x_va);
  compute_bias(ds.nb_samples, ds.y_va, ds.x_va, balanced_error);
}

void LinearClassifier::learn_perceptron(const DataSet &ds, bool balanced_error) {
  for(int i = 0; i < nb_features; i++) weights[i] = 0.0;

  int n_loop_max = 5000;

  for(int i = 0; i < n_loop_max * ds.nb_samples; i++) {
    int ns = i % ds.nb_samples;
    float r = 0;
    for(int f = 0; f < nb_features; f++)
      if(fe_get_bit(ns, ds.x_va[f])) r += weights[f]; else r -= weights[f];
    float correct;
    if(fe_get_bit(ns, ds.y_va)) correct = 1.0; else correct = -1.0;
    if((r < 0 && correct >= 0) || (r >= 0 && correct < 0)) {
      for(int f = 0; f < nb_features; f++)
        if(fe_get_bit(ns, ds.x_va[f])) weights[f] += correct; else weights[f] += -correct;
    }
  }

  compute_bias(ds.nb_samples, ds.y_va, ds.x_va, balanced_error);
}

void FeatureSelector::cmim(const DataSet &ds) {
  fe_selection_cmim(ds.nb_samples, ds.nb_features, ds.x_va, ds.y_va, nb_selected_features, selected_index);
}

void FeatureSelector::mim(const DataSet &ds) {
  fe_selection_mim(ds.nb_samples, ds.nb_features, ds.x_va, ds.y_va, nb_selected_features, selected_index);
}

void FeatureSelector::random(const DataSet &ds) {
  bool used[ds.nb_features];
  for(int i = 0; i < ds.nb_features; i++) used[i] = false;
  int f;
  for(int nf = 0; nf < nb_selected_features; nf++) {
    do { f = int(drand48() * ds.nb_features); } while(used[f]);
    used[f] = true;
    selected_index[nf] = f;
  }
}

FeatureSelector::FeatureSelector(istream &is) {
  is >> nb_selected_features;
  delete[] weights;
  weights = new float[nb_selected_features];
  delete[] selected_index;
  selected_index = new int[nb_selected_features];
  for(int i = 0; i < nb_selected_features; i++) is >> selected_index[i];
}

void FeatureSelector::save(ostream &os) {
  os << nb_selected_features << "\n";
  for(int i = 0; i < nb_selected_features; i++) os << selected_index[i] << ((i < nb_selected_features-1) ? " " : "\n");
}

void LinearClassifier::inner_load(istream &is) {
  is >> nb_features;
  delete[] weights;
  weights = new float[nb_features];
  for(int i = 0; i < nb_features; i++) is >> weights[i];
  is >> bias;
}

void LinearClassifier::save(ostream &os) {
  os << ID_LINEAR << "\n";
  os << nb_features << "\n";
  for(int i = 0; i < nb_features; i++) os << weights[i] << ((i < nb_features-1) ? " " : "\n");
  os << bias << "\n";
}

void LinearClassifier::predict(const DataSet &ds, float *result) {
  for(int s = 0; s < ds.nb_samples; s++) {
    float r = bias;
    for(int nf = 0; nf < nb_features; nf++) if(fe_get_bit(s, ds.x_va[nf])) r += weights[nf];
    result[s] = r;
  }
}

void compute_error_rates(FeatureSelector *selector, Classifier *classifier,
                         const DataSet &testing_set, int &n00, int &n01, int &n10, int &n11, float *result) {

  DataSet reduced_test_set(testing_set, *selector);

  classifier->predict(reduced_test_set, result);

  n00 = 0; n01 = 0; n10 = 0; n11 = 0;
  for(int s = 0; s < testing_set.nb_samples; s++) {
    if(fe_get_bit(s, testing_set.y_va)) {
      if(result[s] >= 0) n11++; else n10++;
    } else {
      if(result[s] >= 0) n01++; else n00++;
    }
  }
}
