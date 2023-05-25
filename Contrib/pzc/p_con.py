# coding=utf-8
# Copyright (c) 2014 Merck KGaA

import csv
import gzip
import json
import math
import optparse
import os
import pickle
import re
import sys
from pickle import Unpickler

import numpy as np
import requests
from pylab import *
from scipy import interp, stats
from sklearn import cross_validation, metrics, preprocessing
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (auc, make_scorer, precision_score, recall_score,
                             roc_curve)

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, SDWriter
from rdkit.ML.Descriptors import MoleculeDescriptors

kappa_template = '''\
%(kind)s Kappa Coefficient
--------------------------------
Kappa %(kappa)6.4f
ASE %(std_kappa)6.4f
%(alpha_ci)s%% Lower Conf Limit %(kappa_low)6.4f
%(alpha_ci)s%% Upper Conf Limit %(kappa_upp)6.4f

Test of H0: %(kind)s Kappa = 0

ASE under H0 %(std_kappa0)6.4f
Z %(z_value)6.4f
One-sided Pr > Z %(pvalue_one_sided)6.4f
Two-sided Pr > |Z| %(pvalue_two_sided)6.4f
'''
'''
Weighted Kappa Coefficient
--------------------------------
Weighted Kappa 0.4701
ASE 0.1457
95% Lower Conf Limit 0.1845
95% Upper Conf Limit 0.7558

Test of H0: Weighted Kappa = 0

ASE under H0 0.1426
Z 3.2971
One-sided Pr > Z 0.0005
Two-sided Pr > |Z| 0.0010
'''


def int_ifclose(x, dec=1, width=4):
  '''helper function for creating result string for int or float

only dec=1 and width=4 is implemented

Parameters
----------
x : int or float
value to format
dec : 1
number of decimals to print if x is not an integer
width : 4
width of string

Returns
-------
xint : int or float
x is converted to int if it is within 1e-14 of an integer
x_string : str
x formatted as string, either '%4d' or '%4.1f'
'''
  xint = int(round(x))
  if np.max(np.abs(xint - x)) < 1e-14:
    return xint, '%4d' % xint
  else:
    return x, '%4.1f' % x


class KappaResults(dict):

  def __init__(self, **kwds):
    self.update(kwds)
    if 'alpha' not in self:
      self['alpha'] = 0.025
      self['alpha_ci'] = int_ifclose(100 - 0.025 * 200)[1]

    self['std_kappa'] = np.sqrt(self['var_kappa'])
    self['std_kappa0'] = np.sqrt(self['var_kappa0'])

    self['z_value'] = self['kappa'] / self['std_kappa0']

    self['pvalue_one_sided'] = stats.norm.sf(self['z_value'])
    self['pvalue_two_sided'] = self['pvalue_one_sided'] * 2

    delta = stats.norm.isf(self['alpha']) * self['std_kappa']
    self['kappa_low'] = self['kappa'] - delta
    self['kappa_upp'] = self['kappa'] + delta

  def __str__(self):
    return kappa_template % self


def cohens_kappa(table, weights=None, return_results=True, wt=None):
  '''Compute Cohen's kappa with variance and equal-zero test

Parameters
----------
table : array_like, 2-Dim
square array with results of two raters, one rater in rows, second
rater in columns
weights : array_like
The interpretation of weights depends on the wt argument.
If both are None, then the simple kappa is computed.
see wt for the case when wt is not None
If weights is two dimensional, then it is directly used as a weight
matrix. For computing the variance of kappa, the maximum of the
weights is assumed to be smaller or equal to one.
TODO: fix conflicting definitions in the 2-Dim case for
wt : None or string
If wt and weights are None, then the simple kappa is computed.
If wt is given, but weights is None, then the weights are set to
be [0, 1, 2, ..., k].
If weights is a one-dimensional array, then it is used to construct
the weight matrix given the following options.

wt in ['linear', 'ca' or None] : use linear weights, Cicchetti-Allison
actual weights are linear in the score "weights" difference
wt in ['quadratic', 'fc'] : use linear weights, Fleiss-Cohen
actual weights are squared in the score "weights" difference
wt = 'toeplitz' : weight matrix is constructed as a toeplitz matrix
from the one dimensional weights.

return_results : bool
If True (default), then an instance of KappaResults is returned.
If False, then only kappa is computed and returned.

Returns
-------
results or kappa
If return_results is True (default), then a results instance with all
statistics is returned
If return_results is False, then only kappa is calculated and returned.

Notes
-----
There are two conflicting definitions of the weight matrix, Wikipedia
versus SAS manual. However, the computation are invariant to rescaling
of the weights matrix, so there is no difference in the results.

Weights for 'linear' and 'quadratic' are interpreted as scores for the
categories, the weights in the computation are based on the pairwise
difference between the scores.
Weights for 'toeplitz' are a interpreted as weighted distance. The distance
only depends on how many levels apart two entries in the table are but
not on the levels themselves.

example:

weights = '0, 1, 2, 3' and wt is either linear or toeplitz means that the
weighting only depends on the simple distance of levels.

weights = '0, 0, 1, 1' and wt = 'linear' means that the first two levels
are zero distance apart and the same for the last two levels. This is
the sampe as forming two aggregated levels by merging the first two and
the last two levels, respectively.

weights = [0, 1, 2, 3] and wt = 'quadratic' is the same as squaring these
weights and using wt = 'toeplitz'.

References
----------
Wikipedia
SAS Manual

'''
  table = np.asarray(table, float)  #avoid integer division
  agree = np.diag(table).sum()
  nobs = table.sum()
  probs = table / nobs
  freqs = probs  #TODO: rename to use freqs instead of probs for observed
  probs_diag = np.diag(probs)
  freq_row = table.sum(1) / nobs
  freq_col = table.sum(0) / nobs
  prob_exp = freq_col * freq_row[:, None]
  assert np.allclose(prob_exp.sum(), 1)
  #print prob_exp.sum()
  agree_exp = np.diag(prob_exp).sum()  #need for kappa_max
  if weights is None and wt is None:
    kind = 'Simple'
    kappa = (agree / nobs - agree_exp) / (1 - agree_exp)

    if return_results:
      #variance
      term_a = probs_diag * (1 - (freq_row + freq_col) * (1 - kappa))**2
      term_a = term_a.sum()
      term_b = probs * (freq_col[:, None] + freq_row)**2
      d_idx = np.arange(table.shape[0])
      term_b[d_idx, d_idx] = 0  #set diagonal to zero
      term_b = (1 - kappa)**2 * term_b.sum()
      term_c = (kappa - agree_exp * (1 - kappa))**2
      var_kappa = (term_a + term_b - term_c) / (1 - agree_exp)**2 / nobs
      #term_c = freq_col * freq_row[:, None] * (freq_col + freq_row[:,None])
      term_c = freq_col * freq_row * (freq_col + freq_row)
      var_kappa0 = (agree_exp + agree_exp**2 - term_c.sum())
      var_kappa0 /= (1 - agree_exp)**2 * nobs

  else:
    if weights is None:
      weights = np.arange(table.shape[0])
    #weights follows the Wikipedia definition, not the SAS, which is 1 -
    kind = 'Weighted'
    weights = np.asarray(weights, float)
    if weights.ndim == 1:
      if wt in ['ca', 'linear', None]:
        weights = np.abs(weights[:, None] - weights) / \
                   (weights[-1] - weights[0])
      elif wt in ['fc', 'quadratic']:
        weights = (weights[:, None] - weights)**2 / \
                   (weights[-1] - weights[0])**2
      elif wt == 'toeplitz':
        #assume toeplitz structure
        from scipy.linalg import toeplitz

        #weights = toeplitz(np.arange(table.shape[0]))
        weights = toeplitz(weights)
      else:
        raise ValueError('wt option is not known')
    else:
      rows, cols = table.shape
      if (table.shape != weights.shape):
        raise ValueError('weights are not square')
    #this is formula from Wikipedia
    kappa = 1 - (weights * table).sum() / nobs / (weights * prob_exp).sum()
    #TODO: add var_kappa for weighted version
    if return_results:
      var_kappa = np.nan
      var_kappa0 = np.nan
      #switch to SAS manual weights, problem if user specifies weights
      #w is negative in some examples,
      #but weights is scale invariant in examples and rough check of source
      w = 1. - weights
      w_row = (freq_col * w).sum(1)
      w_col = (freq_row[:, None] * w).sum(0)
      agree_wexp = (w * freq_col * freq_row[:, None]).sum()
      term_a = freqs * (w - (w_col + w_row[:, None]) * (1 - kappa))**2
      fac = 1. / ((1 - agree_wexp)**2 * nobs)
      var_kappa = term_a.sum() - (kappa - agree_wexp * (1 - kappa))**2
      var_kappa *= fac

      freqse = freq_col * freq_row[:, None]
      var_kappa0 = (freqse * (w - (w_col + w_row[:, None]))**2).sum()
      var_kappa0 -= agree_wexp**2
      var_kappa0 *= fac

  kappa_max = (np.minimum(freq_row, freq_col).sum() - agree_exp) / \
              (1 - agree_exp)

  if return_results:
    res = KappaResults(kind=kind, kappa=kappa, kappa_max=kappa_max, weights=weights,
                       var_kappa=var_kappa, var_kappa0=var_kappa0)
    return res
  else:
    return kappa


def to_table(data, bins=None):
  '''convert raw data with shape (subject, rater) to (rater1, rater2)

    brings data into correct format for cohens_kappa

    Parameters
    ----------
    data : array_like, 2-Dim
        data containing category assignment with subjects in rows and raters
        in columns.
    bins : None, int or tuple of array_like
        If None, then the data is converted to integer categories,
        0,1,2,...,n_cat-1. Because of the relabeling only category levels
        with non-zero counts are included.
        If this is an integer, then the category levels in the data are already
        assumed to be in integers, 0,1,2,...,n_cat-1. In this case, the
        returned array may contain columns with zero count, if no subject
        has been categorized with this level.
        If bins are a tuple of two array_like, then the bins are directly used
        by ``numpy.histogramdd``. This is useful if we want to merge categories.

    Returns
    -------
    arr : nd_array, (n_cat, n_cat)
        Contingency table that contains counts of category level with rater1
        in rows and rater2 in columns.

    Notes
    -----
    no NaN handling, delete rows with missing values

    This works also for more than two raters. In that case the dimension of
    the resulting contingency table is the same as the number of raters
    instead of 2-dimensional.

    '''

  data = np.asarray(data)
  n_rows, n_cols = data.shape
  if bins is None:
    #I could add int conversion (reverse_index) to np.unique
    cat_uni, cat_int = np.unique(data.ravel(), return_inverse=True)
    n_cat = len(cat_uni)
    data_ = cat_int.reshape(data.shape)
    bins_ = np.arange(n_cat + 1) - 0.5
    #alternative implementation with double loop
    #tt = np.asarray([[(x == [i,j]).all(1).sum() for j in cat_uni]
    #                 for i in cat_uni] )
    #other alternative: unique rows and bincount
  elif np.isscalar(bins):
    bins_ = np.arange(bins + 1) - 0.5
    data_ = data
  else:
    bins_ = bins
    data_ = data

  tt = np.histogramdd(data_, (bins_, ) * n_cols)

  return tt[0], bins_


class p_con:
  """Class to create Models to classify Molecules active or inactive
    using threshold for value in training-data"""

  def __init__(self, acc_id=None, proxy={}):
    """Constructor to initialize Object, use proxy if necessary"""
    self.request_data = {"acc_id": acc_id, "proxy": proxy}
    self.acc_id = acc_id
    self.proxy = proxy
    self.model = []
    self.verbous = False

  def __str__(self):
    """String-Representation for Object"""
    self.request_data["cmpd_count"] = len(self.sd_entries)
    retString = ""
    for key in self.request_data.keys():
      retString += "%s: %s\n" % (key, self.request_data[key])
    return retString.rstrip()

  def step_0_get_chembl_data(self):
    """Download Compound-Data for self.acc_id, these are available in self.sd_entries afterwards"""

    def looks_like_number(x):
      """Check for proper Float-Value"""
      try:
        float(x)
        return True
      except ValueError:
        return False

    if self.acc_id.find("CHEMBL") == -1:
      self.target_data = requests.get(
        "https://www.ebi.ac.uk/chemblws/targets/uniprot/{}.json".format(self.acc_id),
        proxies=self.proxy).json()

    else:
      self.target_data = {}
      self.target_data['target'] = {}
      self.target_data['target']['chemblId'] = self.acc_id

    self.chembl_id = self.target_data['target']['chemblId']
    self.request_data["chembl_id"] = self.target_data['target']['chemblId']
    #        print self.target_data
    self.bioactivity_data = requests.get(
      "https://www.ebi.ac.uk/chemblws/targets/{}/bioactivities.json".format(
        self.target_data['target']['chemblId']), proxies=self.proxy).json()

    ic50_skip = 0
    ki_skip = 0
    inhb_skip = 0

    count = 0
    non_homo = 0
    self.dr = {}
    i = 0
    x = len(self.bioactivity_data['bioactivities'])

    for bioactivity in [
        record for record in self.bioactivity_data['bioactivities']
        if looks_like_number(record['value'])
    ]:

      if i % 100 == 0:
        sys.stdout.write('\r' + str(i) + '/' + str(x) + ' >          <\b\b\b\b\b\b\b\b\b\b\b')
      elif (i % 100) % 10 == 0:
        sys.stdout.write('|')
      sys.stdout.flush()
      i += 1
      #            if i > 5000: break
      if bioactivity['organism'] != 'Homo sapiens':
        non_homo += 1
        continue
      if re.search('IC50', bioactivity['bioactivity_type']):
        if bioactivity['units'] != 'nM':
          ic50_skip += 1
          continue
      elif re.search('Ki', bioactivity['bioactivity_type']):
        ki_skip += 1
        continue
      elif re.search('Inhibition', bioactivity['bioactivity_type']):
        inhb_skip += 1
      else:
        continue

      self.cmpd_data = requests.get(
        "https://www.ebi.ac.uk/chemblws/compounds/{}.json".format(
          bioactivity['ingredient_cmpd_chemblid']), proxies=self.proxy).json()

      my_smiles = self.cmpd_data['compound']['smiles']
      bioactivity['Smiles'] = my_smiles
      self.dr[count] = bioactivity
      count += 1

    SDtags = self.dr[0].keys()
    cpd_counter = 0
    self.sd_entries = []
    for x in range(len(self.dr)):
      entry = self.dr[x]
      cpd = Chem.MolFromSmiles(str(entry['Smiles']))
      AllChem.Compute2DCoords(cpd)
      cpd.SetProp("_Name", str(cpd_counter))
      cpd_counter += 1
      for tag in SDtags:
        cpd.SetProp(str(tag), str(entry[tag]))
      self.sd_entries.append(cpd)
    return True

  def step_1_keeplargestfrag(self):
    """remove all smaller Fragments per compound, just keep the largest"""
    result = []

    for cpd in self.sd_entries:
      fragments = Chem.GetMolFrags(cpd, asMols=True)
      list_cpds_fragsize = []
      for frag in fragments:
        list_cpds_fragsize.append(frag.GetNumAtoms())
      largest_frag_index = list_cpds_fragsize.index(max(list_cpds_fragsize))
      largest_frag = fragments[largest_frag_index]
      result.append(largest_frag)

    self.sd_entries = result
    return True

  def step_2_remove_dupl(self):
    """remove duplicates from self.sd_entries"""
    result = []
    all_struct_dict = {}
    for cpd in self.sd_entries:
      Chem.RemoveHs(cpd)
      cansmi = Chem.MolToSmiles(cpd, canonical=True)
      if cansmi not in all_struct_dict.keys():
        all_struct_dict[cansmi] = []
      all_struct_dict[cansmi].append(cpd)

    for entry in all_struct_dict.keys():
      if len(all_struct_dict[entry]) == 1:
        all_struct_dict[entry][0].SetProp('cansmirdkit', entry)
        result.append(all_struct_dict[entry][0])

    self.sd_entries = result
    return True

  def step_3_merge_IC50(self):
    """merge IC50 of duplicates into one compound using mean of all values if:
        min(IC50) => IC50_avg-3*IC50_stddev && max(IC50) <= IC50_avg+3*IC50_stddev && IC50_stddev <= IC50_avg"""
    np_old_settings = np.seterr(invalid='ignore')  #dirty way to ignore warnings from np.std

    def get_mean_IC50(mol_list):
      IC50 = 0
      IC50_avg = 0
      for bla in mol_list:
        try:
          IC50 += float(bla.GetProp("value"))
        except Exception:
          print("no IC50 reported", bla.GetProp("_Name"))
      IC50_avg = IC50 / len(mol_list)
      return IC50_avg

    def get_stddev_IC50(mol_list):
      IC50_list = []
      for mol in mol_list:
        try:
          IC50_list.append(round(float(mol.GetProp("value")), 2))
        except Exception:
          print("no IC50 reported", mol.GetProp("_Name"))
      IC50_stddev = np.std(IC50_list, ddof=1)
      return IC50_stddev, IC50_list

    result = []
    IC50_dict = {}
    for cpd in self.sd_entries:
      if "cansmirdkit" not in cpd.GetPropNames():
        Chem.RemoveHs(cpd)
        cansmi = Chem.MolToSmiles(cpd, canonical=True)
        cpd.SetProp('cansmirdkit', cansmi)
      cansmi = str(cpd.GetProp("cansmirdkit"))
      IC50_dict[cansmi] = {}

    for cpd in self.sd_entries:
      cansmi = str(cpd.GetProp("cansmirdkit"))
      try:
        IC50_dict[cansmi].append(cpd)
      except Exception:
        IC50_dict[cansmi] = [cpd]
    for entry in IC50_dict:
      IC50_avg = str(get_mean_IC50(IC50_dict[entry]))
      IC50_stddev, IC50_list = get_stddev_IC50(IC50_dict[entry])
      IC50_dict[entry][0].SetProp("value_stddev", str(IC50_stddev))
      IC50_dict[entry][0].SetProp("value", IC50_avg)
      minimumvalue = float(IC50_avg) - 3 * float(IC50_stddev)
      maximumvalue = float(IC50_avg) + 3 * float(IC50_stddev)

      if round(IC50_stddev, 1) == 0.0:
        result.append(IC50_dict[entry][0])
      elif IC50_stddev > float(IC50_avg):
        runawaylist = []
        for e in IC50_dict[entry]:
          runawaylist.append(e.GetProp("_Name"))
          print("stddev larger than mean", runawaylist, IC50_list, IC50_avg, IC50_stddev)
      elif np.min(IC50_list) < minimumvalue or np.max(IC50_list) > maximumvalue:
        pass
      else:
        result.append(IC50_dict[entry][0])

    self.sd_entries = result
    np.seterr(over=np_old_settings['over'], divide=np_old_settings['divide'],
              invalid=np_old_settings['invalid'], under=np_old_settings['under'])
    return True

  def step_4_set_TL(self, threshold, ic50_tag="value"):
    """set Property "TL"(TrafficLight) for each compound:
        if ic50_tag (default:"value") > threshold: TL = 0, else 1"""
    result = []
    i, j = 0, 0
    for cpd in self.sd_entries:
      if float(cpd.GetProp(ic50_tag)) > float(threshold):
        cpd.SetProp('TL', '0')
        i += 1
      else:
        cpd.SetProp('TL', '1')
        j += 1
      result.append(cpd)

    self.sd_entries = result
    if self.verbous:
      print("## act: %d, inact: %d" % (j, i))
    return True

  def step_5_remove_descriptors(self):
    """remove list of Properties from each compound (hardcoded)
        which would corrupt process of creating Prediction-Models"""
    sd_tags = [
      'activity__comment', 'alogp', 'assay__chemblid', 'assay__description', 'assay__type',
      'bioactivity__type', 'activity_comment', 'assay_chemblid', 'assay_description', 'assay_type',
      'bioactivity_type', 'cansmirdkit', 'ingredient__cmpd__chemblid', 'ingredient_cmpd_chemblid',
      'knownDrug', 'medChemFriendly', 'molecularFormula', 'name__in__reference',
      'name_in_reference', 'numRo5Violations', 'operator', 'organism', 'parent__cmpd__chemblid',
      'parent_cmpd_chemblid', 'passesRuleOfThree', 'preferredCompoundName', 'reference',
      'rotatableBonds', 'smiles', 'Smiles', 'stdInChiKey', 'synonyms', 'target__chemblid',
      'target_chemblid', 'target__confidence', 'target__name', 'target_confidence', 'target_name',
      'units', 'value_avg', 'value_stddev'
    ] + ['value']
    result = []

    for mol in self.sd_entries:
      properties = mol.GetPropNames()
      for tag in properties:
        if tag in sd_tags:
          mol.ClearProp(tag)
      result.append(mol)

    self.sd_entries = result
    return True

  def step_6_calc_descriptors(self):
    """calculate descriptors for each compound, according to Descriptors._descList"""
    nms = [x[0] for x in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    for i in range(len(self.sd_entries)):
      descrs = calc.CalcDescriptors(self.sd_entries[i])
      for j in range(len(descrs)):
        self.sd_entries[i].SetProp(str(nms[j]), str(descrs[j]))
    return True

  def step_7_train_models(self):
    """train models according to trafficlight using sklearn.ensamble.RandomForestClassifier
        self.model contains up to 10 models afterwards, use save_model_info(type) to create csv or html
        containing data for each model"""
    title_line = [
      "#", "accuracy", "MCC", "precision", "recall", "f1", "auc", "kappa", "prevalence", "bias",
      "pickel-File"
    ]
    self.csv_text = [title_line]

    TL_list = []
    property_list_list = []
    directory = os.getcwd().split("/")[-2:]
    dir_string = ';'.join(directory)
    for cpd in self.sd_entries:
      property_list = []
      property_name_list = []
      prop_name = cpd.GetPropNames()
      for property in prop_name:
        if property not in ['TL', 'value']:
          try:
            f = float(cpd.GetProp(property))
            if math.isnan(f) or math.isinf(f):
              print("invalid: %s" % property)

          except ValueError:
            print("valerror: %s" % property)
            continue
          property_list.append(f)
          property_name_list.append(property)
        elif property == 'TL':
          TL_list.append(int(cpd.GetProp(property)))
        else:
          print(property)
          pass
      property_list_list.append(property_list)
    dataDescrs_array = np.asarray(property_list_list)
    dataActs_array = np.array(TL_list)

    for randomseedcounter in range(1, 11):
      if self.verbous:
        print("################################")
        print("try to calculate seed %d" % randomseedcounter)
      X_train, X_test, y_train, y_test = cross_validation.train_test_split(
        dataDescrs_array, dataActs_array, test_size=.4, random_state=randomseedcounter)
      #            try:
      clf_RF = RandomForestClassifier(n_estimators=100, random_state=randomseedcounter)
      clf_RF = clf_RF.fit(X_train, y_train)

      cv_counter = 5

      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter,
                                                scoring='accuracy')

      accuracy_CV = round(scores.mean(), 3)
      accuracy_std_CV = round(scores.std(), 3)

      calcMCC = make_scorer(metrics.matthews_corrcoef, greater_is_better=True,
                            needs_threshold=False)
      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter,
                                                scoring=calcMCC)

      MCC_CV = round(scores.mean(), 3)
      MCC_std_CV = round(scores.std(), 3)

      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter, scoring='f1')
      scores_rounded = [round(x, 3) for x in scores]
      f1_CV = round(scores.mean(), 3)
      f1_std_CV = round(scores.std(), 3)

      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter,
                                                scoring='precision')
      scores_rounded = [round(x, 3) for x in scores]
      precision_CV = round(scores.mean(), 3)
      precision_std_CV = round(scores.std(), 3)

      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter,
                                                scoring='recall')
      scores_rounded = [round(x, 3) for x in scores]
      recall_CV = round(scores.mean(), 3)
      recall_std_CV = round(scores.std(), 3)

      scores = cross_validation.cross_val_score(clf_RF, X_test, y_test, cv=cv_counter,
                                                scoring='roc_auc')
      scores_rounded = [round(x, 3) for x in scores]
      auc_CV = round(scores.mean(), 3)
      auc_std_CV = round(scores.std(), 3)

      y_predict = clf_RF.predict(X_test)
      conf_matrix = metrics.confusion_matrix(y_test, y_predict)
      #                coh_kappa = cohenskappa.cohens_kappa(conf_matrix)
      coh_kappa = cohens_kappa(conf_matrix)
      kappa = round(coh_kappa['kappa'], 3)
      kappa_stdev = round(coh_kappa['std_kappa'], 3)

      tp = conf_matrix[0][0]
      tn = conf_matrix[1][1]
      fp = conf_matrix[1][0]
      fn = conf_matrix[0][1]
      n = tn + fp
      p = tp + fn
      kappa_prevalence = round(float(abs(tp - tn)) / float(n), 3)
      kappa_bias = round(float(abs(fp - fn)) / float(n), 3)

      if self.verbous:
        print("test:")
        print("\tpos\tneg")
        print("true\t%d\t%d" % (tp, tn))
        print("false\t%d\t%d" % (fp, fn))
        print(conf_matrix)
        print("\ntrain:")
        y_predict2 = clf_RF.predict(X_train)
        conf_matrix2 = metrics.confusion_matrix(y_train, y_predict2)
        tp2 = conf_matrix2[0][0]
        tn2 = conf_matrix2[1][1]
        fp2 = conf_matrix2[1][0]
        fn2 = conf_matrix2[0][1]
        print("\tpos\tneg")
        print("true\t%d\t%d" % (tp2, tn2))
        print("false\t%d\t%d" % (fp2, fn2))
        print(conf_matrix2)

      result_string_cut = [
        randomseedcounter,
        str(accuracy_CV) + "_" + str(accuracy_std_CV),
        str(MCC_CV) + "_" + str(MCC_std_CV),
        str(precision_CV) + "_" + str(precision_std_CV),
        str(recall_CV) + "_" + str(recall_std_CV),
        str(f1_CV) + "_" + str(f1_std_CV),
        str(auc_CV) + "_" + str(auc_std_CV),
        str(kappa) + "_" + str(kappa_stdev), kappa_prevalence, kappa_bias, "model_file.pkl"
      ]

      self.model.append(clf_RF)
      self.csv_text.append(result_string_cut)


#            except Exception as e:
#                print "got %d models" % len(self.model)
#                print e
#                sys.exit(-1)
#                break
    return True if len(self.model) > 0 else False

  def save_model_info(self, outfile, mode="html"):
    """create html- or csv-File for models according to mode (default: "html")"""
    if mode == "csv":
      if not outfile.endswith(".csv"):
        outfile += ".csv"
      csv_file = open(outfile, "wb")
      csv_file_writer = csv.writer(csv_file, delimiter=";", quotechar=' ')
      for line in self.csv_text:
        csv_file_writer.writerow(line)
      csv_file.flush()
      csv_file.close()
    elif mode == "html":
      if not outfile.endswith(".html"):
        outfile += ".html"

      def lines2list(lines):
        return lines

      def list2html(data, act, inact):
        html_head = """<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title></title>
<style type="text/css">
table {
  max-width: 100%;
  background-color: transparent;
}

th {
  text-align: left;
}

.table {
  width: 100%;
  margin-bottom: 20px;
}

.table > thead > tr > th,
.table > tbody > tr > th,
.table > tfoot > tr > th,
.table > thead > tr > td,
.table > tbody > tr > td,
.table > tfoot > tr > td {
  padding: 8px;
  line-height: 1.428571429;
  vertical-align: top;
  border-top: 1px solid #dddddd;
}

.table > thead > tr > th {MSC1013123
  vertical-align: bottom;
  border-bottom: 2px solid #dddddd;
}

.table > caption + thead > tr:first-child > th,
.table > colgroup + thead > tr:first-child > th,
.table > thead:first-child > tr:first-child > th,
.table > caption + thead > tr:first-child > td,
.table > colgroup + thead > tr:first-child > td,
.table > thead:first-child > tr:first-child > td {
  border-top: 0;
}

.table > tbody + tbody {
  border-top: 2px solid #dddddd;
}

.table .table {
  background-color: #ffffff;
}

.table-condensed > thead > tr > th,
.table-condensed > tbody > tr > th,
.table-condensed > tfoot > tr > th,
.table-condensed > thead > tr > td,
.table-condensed > tbody > tr > td,
.table-condensed > tfoot > tr > td {
  padding: 5px;
}

.table-bordered {
  border: 1px solid #dddddd;
}

.table-bordered > thead > tr > th,
.table-bordered > tbody > tr > th,
.table-bordered > tfoot > tr > th,
.table-bordered > thead > tr > td,
.table-bordered > tbody > tr > td,
.table-bordered > tfoot > tr > td {
  border: 1px solid #dddddd;
}

.table-bordered > thead > tr > th,
.table-bordered > thead > tr > td {
  border-bottom-width: 2px;
}

.table-striped > tbody > tr:nth-child(odd) > td,
.table-striped > tbody > tr:nth-child(odd) > th {
  background-color: #f9f9f9;
}

.table-hover > tbody > tr:hover > td,
.table-hover > tbody > tr:hover > th {
  background-color: #f5f5f5;
}

table col[class*="col-"] {
  position: static;
  display: table-column;
  float: none;
}

table td[class*="col-"],
table th[class*="col-"] {
  display: table-cell;
  float: none;
}

.table > thead > tr > .active,
.table > tbody > tr > .active,
.table > tfoot > tr > .active,
.table > thead > .active > td,
.table > tbody > .active > td,
.table > tfoot > .active > td,
.table > thead > .active > th,
.table > tbody > .active > th,
.table > tfoot > .active > th {
  background-color: #f5f5f5;
}

.table-hover > tbody > tr > .active:hover,
.table-hover > tbody > .active:hover > td,
.table-hover > tbody > .active:hover > th {
  background-color: #e8e8e8;
}

.table > thead > tr > .success,
.table > tbody > tr > .success,
.table > tfoot > tr > .success,
.table > thead > .success > td,
.table > tbody > .success > td,
.table > tfoot > .success > td,
.table > thead > .success > th,
.table > tbody > .success > th,
.table > tfoot > .success > th {
  background-color: #dff0d8;
}

.table-hover > tbody > tr > .success:hover,
.table-hover > tbody > .success:hover > td,
.table-hover > tbody > .success:hover > th {
  background-color: #d0e9c6;
}

.table > thead > tr > .danger,
.table > tbody > tr > .danger,
.table > tfoot > tr > .danger,
.table > thead > .danger > td,
.table > tbody > .danger > td,
.table > tfoot > .danger > td,
.table > thead > .danger > th,
.table > tbody > .danger > th,
.table > tfoot > .danger > th {
  background-color: #f2dede;
}

.table-hover > tbody > tr > .danger:hover,
.table-hover > tbody > .danger:hover > td,
.table-hover > tbody > .danger:hover > th {
  background-color: #ebcccc;
}

.table > thead > tr > .warning,
.table > tbody > tr > .warning,
.table > tfoot > tr > .warning,
.table > thead > .warning > td,
.table > tbody > .warning > td,
.table > tfoot > .warning > td,
.table > thead > .warning > th,
.table > tbody > .warning > th,
.table > tfoot > .warning > th {
  background-color: #fcf8e3;
}

.table-hover > tbody > tr > .warning:hover,
.table-hover > tbody > .warning:hover > td,
.table-hover > tbody > .warning:hover > th {
  background-color: #faf2cc;
}

@media (max-width: 767px) {
  .table-responsive {
    width: 100%;
    margin-bottom: 15px;
    overflow-x: scroll;
    overflow-y: hidden;
    border: 1px solid #dddddd;
    -ms-overflow-style: -ms-autohiding-scrollbar;
    -webkit-overflow-scrolling: touch;
  }
  .table-responsive > .table {
    margin-bottom: 0;
  }
  .table-responsive > .table > thead > tr > th,
  .table-responsive > .table > tbody > tr > th,
  .table-responsive > .table > tfoot > tr > th,
  .table-responsive > .table > thead > tr > td,
  .table-responsive > .table > tbody > tr > td,
  .table-responsive > .table > tfoot > tr > td {
    white-space: nowrap;
  }
  .table-responsive > .table-bordered {
    border: 0;
  }
  .table-responsive > .table-bordered > thead > tr > th:first-child,
  .table-responsive > .table-bordered > tbody > tr > th:first-child,
  .table-responsive > .table-bordered > tfoot > tr > th:first-child,
  .table-responsive > .table-bordered > thead > tr > td:first-child,
  .table-responsive > .table-bordered > tbody > tr > td:first-child,
  .table-responsive > .table-bordered > tfoot > tr > td:first-child {
    border-left: 0;
  }
  .table-responsive > .table-bordered > thead > tr > th:last-child,
  .table-responsive > .table-bordered > tbody > tr > th:last-child,
  .table-responsive > .table-bordered > tfoot > tr > th:last-child,
  .table-responsive > .table-bordered > thead > tr > td:last-child,
  .table-responsive > .table-bordered > tbody > tr > td:last-child,
  .table-responsive > .table-bordered > tfoot > tr > td:last-child {
    border-right: 0;
  }
  .table-responsive > .table-bordered > tbody > tr:last-child > th,
  .table-responsive > .table-bordered > tfoot > tr:last-child > th,
  .table-responsive > .table-bordered > tbody > tr:last-child > td,
  .table-responsive > .table-bordered > tfoot > tr:last-child > td {
    border-bottom: 0;
  }
}
</style>
</head>
<body>
<p style="padding-left:10px;padding-top:10px;font-size:200&#37;">Data for Models</p>
<p style="padding-left:10px;padding-right:10px;">"""

        html_topPlot_start = """<table style="vertical-align:top; background-color=#CCCCCC">
<tr align="left" valign="top"><td><img src="pieplot.png"></td><td><H3>Distribution</H3><font color="#00C000">active %d</font><br><font color="#FF0000">inactive %d</td><td>"""

        html_topPlot_bottom = """</td></tr></table>"""

        html_tableStart = """<table class="table table-bordered table-condensed">
<thead>
<tr>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
</tr>
</thead>
<tbody>"""

        html_tElements = """
<tr bgcolor = "%s">
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td><a href="%s">model.pkl</a></td>
</tr>"""

        html_bottomPlot = """</tbody>
</table>
<img src="barplot.png"><br>"""

        html_foot = """
</p>
</body>
</html>"""

        html_kappa_table_head = """<table class="table table-bordered table-condensed">
<thead>
<tr>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
<th>%s</th>
</tr>
</thead>
<tbody>"""

        html_kappa_table_element = """<tr bgcolor = "%s">
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td>%s</td>
<td><a href="%s">model.pkl</a></td>
</tr>"""

        html_kappa_table_bottom = """</tbody>
</table>
<img src="barplot.png"><br>"""

        best, worst = findBestWorst(data)
        html = []
        html.append(html_head)
        html.append(html_topPlot_start % (act, inact))
        html.append(html_topPlot_bottom)
        html.append(html_tableStart % tuple(data[0]))
        i = 0
        for l in data[1:len(data)]:
          l_replaced = []
          for elem in l:
            elem_string = str(elem)
            if elem_string.find("pkl") == -1:
              l_replaced.append(elem_string.replace("_", "±"))
            else:
              l_replaced.append(elem_string)

          c = ""
          if i == best:
            c = "#9CC089"
          if i == worst:
            c = "#FF3333"

          html.append(html_tElements % tuple([c] + l_replaced))
          i += 1
        html.append(html_bottomPlot)
        html.append(html_foot)
        createBarPlot(data)
        return html

      def writeHtml(html, outf):
        outf_h = open(outf, 'w')
        for block in html:
          outf_h.write(block)
        outf_h.flush()
        outf_h.close()
        return

      def findBestWorst(data):
        auc = [float(x[6].split("_")[0]) for x in data[1:]]
        max_index, min_index = auc.index(max(auc)), auc.index(min(auc))
        return (max_index, min_index)

      def createPiePlot(cpds):

        def getActInact(cpds):
          act, inact = 0, 0
          for cpd in cpds:
            if int(cpd.GetProp('TL')) == 0:
              inact += 1
            else:
              act += 1
          return act, inact

        act_count, inact_count = getActInact(cpds)
        print("act/inact from TL's %d/%d" % (act_count, inact_count))
        fig = plt.figure(figsize=(2, 2))
        pie = plt.pie([inact_count, act_count], colors=('r', 'g'))
        fig.savefig("pieplot.png", transparent=True)
        return act_count, inact_count

      def createBarPlot(data):

        def getLists(data, col):
          accList = []
          errList = []
          for x in data[1:]:
            if x[col].find("_") == -1:
              continue
            if x[col].find(".pkl") != -1:
              continue
            spl = x[col].split("_")
            accList.append(float(spl[0]))
            errList.append(float(spl[1]))
          return accList, errList

        def plotLists(cnt):
          result = []
          clr = ['#DD1E2F', '#EBB035', '#06A2CB', '#218559', '#D0C6B1', '#192823', '#DDAACC']
          #                    print ticks, list,errList,width
          #                    print ticks
          for i in range(1, cnt):
            list, errList = getLists(data, i)
            #                        print i,cnt,list,errList
            result.append(ax.bar(ticks + width * i, list, width, color=clr[i - 1], yerr=errList))
          return result

        fig, ax = plt.subplots()
        fig.set_size_inches(15, 6)
        ticks = np.arange(0.0, 12.0, 1.2)
        if len(self.model) == 1:
          ticks = np.arange(0.0, 1.0, 1.5)
        width = 0.15
        plots = plotLists(8)
        ax.set_xticks(ticks + 0.75)
        ax.set_xticklabels([str(x) for x in range(1, 11, 1)])
        ax.set_ylabel("Accuracy")
        ax.set_xlabel("# model")
        ax.set_xlim(-0.3, 14)
        ax.set_ylim(-0.1, 1.2)
        ax.legend(tuple(plots), [x for x in data[0][1:8]], 'upper right')
        best, worst = findBestWorst(data)
        if len(self.model) > 1:
          ax.annotate("best", xy=(ticks[best], 0.85), xytext=(ticks[best] + 0.25, 1.1),
                      color="green")
          ax.annotate("worst", xy=(ticks[worst], 0.85), xytext=(ticks[worst] + 0.25, 1.10),
                      color="red")
        fig.savefig("barplot.png", transparent=True)
        return

      act, inact = createPiePlot(self.sd_entries)
      lines = self.csv_text
      data = lines2list(lines)
      html = list2html(data, act, inact)
      writeHtml(html, outfile)
    return True

  def load_mols(self, sd_file):
    """load SD-File from .sdf, .sdf.gz or .sd.gz"""
    if sd_file.endswith(".sdf.gz") or sd_file.endswith(".sd.gz"):
      SDFile = Chem.ForwardSDMolSupplier(gzip.open(sd_file))
    else:
      SDFile = Chem.SDMolSupplier(sd_file)
    self.sd_entries = [mol for mol in SDFile]
    return True

  def save_mols(self, outfile, gzip=True):
    """create SD-File of current molecules in self.sd_entries"""
    sdw = Chem.SDWriter(outfile + ".tmp")
    for mol in self.sd_entries:
      sdw.write(mol)
    sdw.flush()
    sdw.close()
    if not gzip:
      os.rename(outfile + ".tmp", outfile)
      return
    f_in = open(outfile + ".tmp", 'rb')
    f_out = gzip.open(outfile, 'wb')
    f_out.writelines(f_in)
    f_out.flush()
    f_out.close()
    f_in.close()
    os.remove(outfile + ".tmp")
    return

  def save_model(self, outfile, model_number=0):
    """save Model to file using pickle.dump"""
    pickle.dump(self.model[model_number], file(outfile, "wb+"))
    return

  def load_models(self, model_files):
    """load model or list of models into self.model"""
    if type(model_files) == str:
      model_files = [model_files]
    i = 0
    for mod_file in model_files:
      model = open(mod_file, 'r')
      unPickled = Unpickler(model)
      clf_RF = unPickled.load()
      self.model.append(clf_RF)
      model.close()
      i += 1
    return i

  def predict(self, model_number):
    """try to predict activity of compounds using giving model-Number"""
    if len(self.model) <= model_number:
      sys.stderr.write("\nModel-Number %d doesn't exist, there are just %d Models\n" %
                       (model_number, len(self.model)))
      sys.exit(-1)
    descriptors = []
    active, inactive = 0, 0

    for D in Descriptors._descList:
      descriptors.append(D[0])
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)

    clf_RF = self.model[model_number]

    for sample in self.sd_entries:
      use = False
      try:
        pattern = calculator.CalcDescriptors(sample)
        use = True
      except e:
        sys.stderr.write("Error computing descriptors for %s, skip" % sample)

      if use:
        dataDescrs_array = np.asarray(pattern)
        y_predict = int(clf_RF.predict(dataDescrs_array)[0])
        if y_predict == 0:
          inactive += 1
        if y_predict == 1:
          active += 1
        sample.SetProp("TL_prediction", str(y_predict))
    return (active, inactive)

if __name__ == "__main__":

  def step_error(step):
    sys.stderr.write("Error in Step: %s" % step)

  usage = "usage: python master.py [--accession=<Acc_ID>] [--sdf=<sdf-File>] --dupl/--uniq [--rof] [--combine=<file1>,<file2>] [--IC50=<IC50_tag>] [--cutoff=<value>] [--remove_descr=<txt_file>] [--proxy=<https://user:pass@proxy.de:portnumber] [--verbous] [--check_models=<model.pkl>]"
  parser = optparse.OptionParser(usage=usage)
  parser.add_option(
    '--accession', action='store', type='string', dest='accession',
    help="Accession ID of Protein (hint: P43088 is Vitamin_D_Receptor with ~200 compounds)",
    default='')
  parser.add_option('--rof', action='store_true', dest='onefile', help='remove obsolete Files',
                    default=False)
  parser.add_option('--dupl', action='store_true', dest='dupl', help='use only duplicates',
                    default=False)
  parser.add_option('--uniq', action='store_true', dest='uniq', help='use only uniques',
                    default=False)
  parser.add_option('--combine', action='store', type='string', dest='combine',
                    help='Combine 2 SDF/SDF.GZ Files', default='')
  parser.add_option('--IC50', action='store', type='string', dest='SD_tag',
                    help='name of IC50 field, default is \'value\'', default='value')
  parser.add_option('--cutoff', action='store', type='int', dest='cutoff',
                    help='cutoff-value for hERG-trafficlight, default is \'5000\'', default=5000)
  parser.add_option('--remove_descr', action='store', type='string', dest='remove_descr',
                    help='file with SDtags2remove, line-wise default:<internal list>', default='')
  parser.add_option('--proxy', action='store', type='string', dest='proxy', help='Use this Proxy',
                    default='')
  parser.add_option('--sdf', action='store', type='string', dest='sdf', help='load this SDF-File',
                    default='')
  parser.add_option('--verbous', action='store_true', dest='verbous', help='verbous', default=False)
  parser.add_option('--check_models', action='store', type='string', dest='modelfile',
                    help='check compounds with this model', default='')

  (options, args) = parser.parse_args()
  combineItems = options.combine.split(',')

  if len(combineItems) == 1 and len(combineItems[0]) > 0:
    print('need 2 files to combine')
    print(usage)
    sys.exit(-1)
  elif len(combineItems) == 2 and len(combineItems[0]) > 0 and len(combineItems[1]) > 0:
    cur_file = _04.combine(combineItems[0], combineItems[1])
    print("File: %s" % cur_file)
    sys.exit(0)

  code = options.accession.split(':')
  if len(code) == 1:
    accession = code[0]
  else:
    accession = code[1]

  if options.accession == '' and options.sdf == '':
    print("please offer Accession-Number or SDF-File")
    print("-h for help")
    sys.exit(-1)

  if not options.dupl and not options.uniq:
    print("Please select uniq or dupl -h for help")
    print("-h for help")
    sys.exit(-1)

  pco = p_con(accession, proxy=options.proxy)
  pco.verbous = options.verbous

  if options.sdf != '':
    print("load sdf from File: %s" % options.sdf)
    result = pco.load_mols(options.sdf)
    if not result:
      step_error("load SDF-File")
      sys.exit(-1)
  else:
    print("gather Data for Accession-ID \'%s\'" % accession)
    result = pco.step_0_get_chembl_data()
    if not result:
      step_error("download ChEMBL-Data")
      sys.exit(-1)

  result = pco.step_1_keeplargestfrag()
  if not result:
    step_error("keep largest Fragment")
    sys.exit(-1)

  if options.uniq:
    result = pco.step_2_remove_dupl()
    if not result:
      step_error("remove duplicates")
      sys.exit(-1)

  result = pco.step_3_merge_IC50()
  if not result:
    step_error("merge IC50-Values for same Smiles")
    sys.exit(-1)

  if options.modelfile != '':
    result = pco.load_models(options.modelfile.split(","))
    if not result:
      step_error("Load Model-Files")
      sys.exit(-1)
    print("\n#Model\tActive\tInactive")
    for i in range(len(pco.model)):
      act, inact = pco.predict(i)
      print("%d\t%d\t%d" % (i, act, inact))
    sys.exit(0)

  result = pco.step_4_set_TL(options.cutoff)
  if not result:
    step_error("set Trafficlight for cutoff")
    sys.exit(-1)

  result = pco.step_5_remove_descriptors()
  if not result:
    step_error("remove descriptors")
    sys.exit(-1)

  result = pco.step_6_calc_descriptors()
  if not result:
    step_error("calculate Descriptors")
    sys.exit(-1)

  result = pco.step_7_train_models()
  if not result:
    step_error("Training of Models")
    sys.exit(-1)

  pco.save_model_info("model_info.csv", mode="csv")
  pco.save_model_info("model_info.html", mode="html")

  for i in range(len(pco.model)):
    filename = "%s_%dnm_model_%d.pkl" % (accession, options.cutoff, i)
    pco.save_model(filename, i)
    print("Model %d saved into File: %s" % (i, filename))

  for i in range(len(pco.model)):
    act, inact = pco.predict(i)
    print("Model %d active: %d\tinactive: %d" % (i, act, inact))
