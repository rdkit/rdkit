
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

// $Id: cmim.cc,v 1.2 2005/03/03 20:16:15 fleuret Exp $

// This software was developped on GNU/Linux systems with many GPL
// tools including emacs, gcc, gdb, and bash (see http://www.fsf.org).

/*

To test

./cmim --feature-selection cmim --classifier bayesian --error ber --train ./train.dat ./classifier.nb 100
./cmim --test ./test.dat ./classifier.nb ./result.dat

*/

using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "classifier.h"

#define BUFFER_SIZE 256

FeatureSelector *selector;
Classifier *classifier;
char classifier_type[BUFFER_SIZE] = "bayesian";
char feature_selection_type[BUFFER_SIZE] = "cmim";
float reg_param = 0.0;
bool verbose = true;
bool balanced_error = false;
int nb_selected_features = 100;

void check_opt(int argc, char **argv, int n_opt, int n, char *help) {
  if(n_opt+n >= argc) {
    cerr << "Missing argument for " << argv[n_opt] << ".\n";
    cerr << "Expecting " << help << ".\n";
    exit(1);
  }
}

void train(const DataSet &training_set) {
  timeval tv_start, tv_end;
  fe_init_tables();

  if(verbose) {
    cout << "Selecting features with " << feature_selection_type;
    cout.flush();
    gettimeofday(&tv_start, 0);
  }

  cout.flush();

  selector = new FeatureSelector(nb_selected_features);

  if(strcmp(feature_selection_type, "cmim") == 0) selector->cmim(training_set);
  else if(strcmp(feature_selection_type, "mim") == 0) selector->mim(training_set);
  else if(strcmp(feature_selection_type, "random") == 0) selector->random(training_set);
  else {
    cerr << "Unknown feature selection type " << feature_selection_type << "\n";
    exit(1);
  }

  if(verbose) {
    gettimeofday(&tv_end, 0);
    cout << " (" <<
      (float(tv_end.tv_sec * 1000) + float(tv_end.tv_usec)/1000) -
      (float(tv_start.tv_sec * 1000) + float(tv_start.tv_usec)/1000) << "ms).\n";

    //for(int i=0;i<nb_selected_features;++i){
    //  cout << selector->selected_index[i] << ", ";
    //}
    //cout << "\n";

    gettimeofday(&tv_start, 0);
    cout << "Learning with " << classifier_type;
    cout.flush();
  }

  cout.flush();

  DataSet reduced_training_set(training_set, *selector);

  if(strcmp(classifier_type, "bayesian") == 0) {
    LinearClassifier *tmp = new LinearClassifier(nb_selected_features);
    tmp->learn_bayesian(reduced_training_set, balanced_error);
    classifier = tmp;
  }

  else if(strcmp(classifier_type, "perceptron") == 0) {
    LinearClassifier *tmp = new LinearClassifier(nb_selected_features);
    tmp->learn_perceptron(reduced_training_set, balanced_error);
    classifier = tmp;
  }

  else {
    cerr << "Unknown learning method type " << classifier_type << "\n";
    exit(1);
  }

  if(verbose) {
    gettimeofday(&tv_end, 0);
    cout << " (" <<
      (float(tv_end.tv_sec * 1000) + float(tv_end.tv_usec)/1000) -
      (float(tv_start.tv_sec * 1000) + float(tv_start.tv_usec)/1000) << "ms).\n";
  }

  cout.flush();
}

int main(int argc, char **argv) {
  bool arg_error = false;
  bool delim=false;

  int i = 1;
  while(i < argc && !arg_error) {

    //////////////////////////////////////////////////////////////////////
    // Parameters ////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    if(strcmp(argv[i], "--silent") == 0) {
      verbose = false;
      i++;
    }

    else if(strcmp(argv[i], "--feature-selection") == 0) {
      check_opt(argc, argv, i, 1, "<random|mim|cmim>");
      strncpy(feature_selection_type, argv[i+1], BUFFER_SIZE);
      i += 2;
    }

    else if(strcmp(argv[i], "--classifier") == 0) {
      check_opt(argc, argv, i, 1, "<bayesian|perceptron>");
      strncpy(classifier_type, argv[i+1], BUFFER_SIZE);
      i += 2;
    }

    else if(strcmp(argv[i], "--error") == 0) {
      check_opt(argc, argv, i, 1, "<standard|ber>");
      if(strcmp(argv[i+1], "standard") == 0) balanced_error = false;
      else if(strcmp(argv[i+1], "ber") == 0) balanced_error = true;
      else {
        cerr << "Unknown  error type " << argv[i+1] << "!\n";
        exit(1);
      }
      i += 2;
    }

    else if(strcmp(argv[i], "--nb-features") == 0) {
      check_opt(argc, argv, i, 1, "<int: nb features>");
      nb_selected_features = atoi(argv[i+1]);
      
      if(nb_selected_features <= 0) {
        cerr << "Unconsistent number of selected features (" << nb_selected_features << ").\n";
        exit(1);
      }
      i += 2;
    }

    //////////////////////////////////////////////////////////////////////
    // Training //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    else if(strcmp(argv[i], "--cross-validation") == 0) {
      check_opt(argc, argv, i, 3, "<file: data set> <int: nb test samples> <int: nb loops>");
      if(verbose) {
        cout << "Loading data.\n";
        cout.flush();
      }

      ifstream complete_data(argv[i+1]);
      if(complete_data.fail()) {
        cerr << "Can not open " << argv[i+1] << " for reading!\n";
        exit(1);
      }

      int nb_for_test = atoi(argv[i+2]);
      if(nb_for_test <= 0) {
        cerr << "Unconsistent number of samples for test (" << nb_selected_features << ").\n";
        exit(1);
      }

      int nb_cv_loops = atoi(argv[i+3]);
      if(nb_cv_loops <= 0) {
        cerr << "Unconsistent number of cross-validation loops (" << nb_cv_loops << ").\n";
        exit(1);
      }

      DataSet complete_set(complete_data,delim);

      int n00_test = 0, n01_test = 0, n10_test = 0, n11_test = 0;
      int n00_train = 0, n01_train = 0, n10_train = 0, n11_train = 0;

      for(int ncv = 0; ncv < nb_cv_loops; ncv++) {
        bool for_test[complete_set.nb_samples];
	cerr << " Round " << ncv << " of " << nb_cv_loops << " " << complete_set.nb_samples << " points \n";cerr.flush();

	cerr << " 1 " << "\n";cerr.flush();
        for(int s = 0; s < complete_set.nb_samples; s++) for_test[s] = false;
	cerr << " 2 " << "\n";cerr.flush();
        for(int i = 0; i < nb_for_test; i++) {
          int s;
          do {
	    double tmpV=drand48();
            s = int(tmpV * complete_set.nb_samples);
          } while(for_test[s]);
          for_test[s] = true;
        }

	cerr << " 3 " << "\n";cerr.flush();
        DataSet testing_set(complete_set, for_test);
	cerr << " 4 " << "\n";cerr.flush();
        for(int s = 0; s < complete_set.nb_samples; s++) for_test[s] = !for_test[s];
	cerr << " 5 " << "\n";cerr.flush();
        DataSet training_set(complete_set, for_test);

	cerr << " training " << "\n";cerr.flush();
        train(training_set);

        int n00, n01, n10, n11;

	cerr << " error rates 1 " << "\n";cerr.flush();
        {
          float result[training_set.nb_samples];
          compute_error_rates(selector, classifier, training_set, n00, n01, n10, n11, result);
          n00_train += n00; n01_train += n01; n10_train += n10; n11_train += n11;
        }

	cerr << " error rates 2 " << "\n";cerr.flush();
        {
          float result[testing_set.nb_samples];
          compute_error_rates(selector, classifier, testing_set, n00, n01, n10, n11, result);
          n00_test += n00; n01_test += n01; n10_test += n10; n11_test += n11;
        }

        delete classifier;
        delete selector;
      }

      if(balanced_error) {
        cout << "BER [" << nb_cv_loops << " loops] "
                  << " train " << 0.5 * (float(n01_train)/float(n00_train + n01_train) + float(n10_train)/float(n10_train + n11_train))
                  << " test " << 0.5 * (float(n01_test)/float(n00_test + n01_test) + float(n10_test)/float(n10_test + n11_test)) << "\n";
      } else {
        cout << "Error [" << nb_cv_loops << " loops] "
                  << " train " << float(n01_train + n10_train)/float(n00_train + n01_train + n10_train + n11_train)
                  << " test " << float(n01_test + n10_test)/float(n00_test + n01_test + n10_test + n11_test) << "\n";
      }

      i += 4;
    }

    //////////////////////////////////////////////////////////////////////

    else if(strcmp(argv[i], "--train") == 0) {
      check_opt(argc, argv, i, 2, "<file: data set> <file: classifier>");

      if(verbose) {
        cout << "Loading data.\n";
        cout.flush();
      }

      ifstream training_data(argv[i+1]);
      if(training_data.fail()) {
        cerr << "Can not open " << argv[i+1] << " for reading!\n";
        exit(1);
      }

      DataSet training_set(training_data,delim);

      //////////////////////////////////////////////////////////////////////
      // Learning with CMIM + naive Bayesian ///////////////////////////////
      //////////////////////////////////////////////////////////////////////

      train(training_set);

      //////////////////////////////////////////////////////////////////////
      // Finishing and saving //////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////

      if(verbose) cout << "Saving the classifier in [" << argv[i+2] << "].\n";
      ofstream classifier_out(argv[i+2]);
      if(classifier_out.fail()) {
        cerr << "Can not open " << argv[i+2] << " for writing!\n";
        exit(1);
      }

      selector->save(classifier_out);
      classifier->save(classifier_out);

      delete classifier;
      delete selector;

      i += 3;
    }

    //////////////////////////////////////////////////////////////////////
    // Test //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    else if(strcmp(argv[i], "--test") == 0) {
      check_opt(argc, argv, i, 3, "<file: classifier> <file: data set> <file: result>");

      // Load the classifier

      if(verbose) cout << "Loading the classifier from [" << argv[i+1] << "].\n";

      ifstream classifier_in(argv[i+1]);
      if(classifier_in.fail()) {
        cerr << "Can not open " << argv[i+1] << " for reading!\n";
        exit(1);
      }

      selector = new FeatureSelector(classifier_in);
      classifier = Classifier::load(classifier_in);

      // Load the testing data

      ifstream testing_data(argv[i+2]);
      if(testing_data.fail()) {
        cerr << "Can not open " << argv[i+2] << " for reading!\n";
        exit(1);
      }

      ofstream result_out(argv[i+3]);
      if(result_out.fail()) {
        cerr << "Can not open " << argv[i+3] << " for writing!\n";
        exit(1);
      }

      DataSet testing_set(testing_data,delim);

      // Compute the predicted responses

      int n00, n01, n10, n11;
      float result[testing_set.nb_samples];
      compute_error_rates(selector, classifier, testing_set, n00, n01, n10, n11, result);

      for(int s = 0; s < testing_set.nb_samples; s++)
        result_out << result[s] << "\n";

      cout << "ERROR " << float(n01 + n10)/float(n00 + n01 + n10 + n11) << "\n";
      cout << "BER   " << 0.5 * (float(n01)/float(n00 + n01) + float(n10)/float(n10 + n11)) << "\n";
      cout << "FN    " << float(n10)/float(n10+n11) << "\n";
      cout << "FP    " << float(n01)/float(n01+n00) << "\n";
      cout << "real_0_predicted_0 " << n00 << "\n";
      cout << "real_0_predicted_1 " << n01 << "\n";
      cout << "real_1_predicted_0 " << n10 << "\n";
      cout << "real_1_predicted_1 " << n11 << "\n";

      delete classifier;
      delete selector;

      i += 4;
    }

    else arg_error = true;
  }

  if(arg_error) {
    cerr << "Conditional Mutual Information Maximization\n";
    cerr << "Written by François Fleuret (c) EPFL 2004\n";
    cerr << "Comments and bug reports to <francois.fleuret@epfl.ch>\n";
    cerr << "\n";
    cerr << "Usage: " << argv[0] << "\n";
    cerr << "--silent\n";
    cerr << "--feature-selection <random|mim|cmim>\n";
    cerr << "--classifier <bayesian|perceptron>\n";
    cerr << "--error <standard|ber>\n";
    cerr << "--nb-features <int: nb of features>\n";
    cerr << "--cross-validation <file: data set> <int: nb test samples> <int: nb loops>\n";
    cerr << "--train <file: data set> <file: classifier>\n";
    cerr << "--test <file: classifier> <file: data set> <file: result>\n";
    exit(1);
  }
}
