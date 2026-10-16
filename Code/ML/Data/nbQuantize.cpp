//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>
#include <typeinfo>
#include <vector>

#include <RDGeneral/Invariant.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

template <typename T>
std::vector<T> sequenceToVector(const nb::object &sequence) {
  auto size = nb::len(sequence);
  std::vector<T> result;
  result.reserve(size);
  for (decltype(size) i = 0; i < size; ++i) {
    try {
      result.push_back(nb::cast<T>(sequence[i]));
    } catch (const std::bad_cast &) {
      throw nb::type_error("expected a one-dimensional sequence of numbers");
    }
  }
  return result;
}

}  // namespace

long int *GenVarTable(double *vals, int nVals, long int *cuts, int nCuts,
                      long int *starts, long int *results, int nPossibleRes,
                      long int *varTable) {
  RDUNUSED_PARAM(vals);
  int nBins = nCuts + 1;
  int idx = 0;
  int iTab;

  memset(varTable, 0, nBins * nPossibleRes * sizeof(long int));
  for (int i = 0; i < nCuts; i++) {
    int cut = cuts[i];
    iTab = i * nPossibleRes;
    while (idx < starts[cut]) {
      varTable[iTab + results[idx]] += 1;
      idx++;
    }
  }
  iTab = nCuts * nPossibleRes;
  while (idx < nVals) {
    varTable[iTab + results[idx]] += 1;
    idx++;
  }
  return varTable;
}

double RecurseHelper(double *vals, int nVals, long int *cuts, int nCuts,
                     int which, long int *starts, int nStarts,
                     long int *results, int nPossibleRes) {
  PRECONDITION(vals, "bad vals pointer");

  double maxGain = -1e6, gainHere;
  int highestCutHere = nStarts - nCuts + which;
  int nBounds = nCuts;

  std::vector<long int> varTable((nCuts + 1) * nPossibleRes);
  std::vector<long int> bestCuts(nCuts);
  std::vector<long int> tCuts(nCuts);
  GenVarTable(vals, nVals, cuts, nCuts, starts, results, nPossibleRes,
              varTable.data());
  while (cuts[which] <= highestCutHere) {
    gainHere =
        RDInfoTheory::InfoEntropyGain(varTable.data(), nCuts + 1, nPossibleRes);
    if (gainHere > maxGain) {
      maxGain = gainHere;
      memcpy(bestCuts.data(), cuts, nCuts * sizeof(long int));
    }

    if (which < nBounds - 1) {
      memcpy(tCuts.data(), cuts, nCuts * sizeof(long int));
      gainHere = RecurseHelper(vals, nVals, tCuts.data(), nCuts, which + 1,
                               starts, nStarts, results, nPossibleRes);
      if (gainHere > maxGain) {
        maxGain = gainHere;
        memcpy(bestCuts.data(), tCuts.data(), nCuts * sizeof(long int));
      }
    }

    int oldCut = cuts[which];
    cuts[which] += 1;
    int top, bot;
    bot = starts[oldCut];
    if (oldCut + 1 < nStarts) {
      top = starts[oldCut + 1];
    } else {
      top = starts[nStarts - 1];
    }
    for (int i = bot; i < top; i++) {
      int v = results[i];
      varTable[which * nPossibleRes + v] += 1;
      varTable[(which + 1) * nPossibleRes + v] -= 1;
    }
    for (int i = which + 1; i < nBounds; i++) {
      if (cuts[i] == cuts[i - 1]) {
        cuts[i] += 1;
      }
    }
  }
  memcpy(cuts, bestCuts.data(), nCuts * sizeof(long int));
  return maxGain;
}

static nb::tuple cQuantize_RecurseOnBounds(nb::object vals, nb::list pyCuts,
                                           int which, nb::list pyStarts,
                                           nb::object results,
                                           int nPossibleRes) {
  auto values = sequenceToVector<double>(vals);
  auto resultCodes = sequenceToVector<long int>(results);
  auto cuts = sequenceToVector<long int>(pyCuts);
  auto starts = sequenceToVector<long int>(pyStarts);
  if (resultCodes.size() < values.size()) {
    throw nb::value_error("results array is shorter than values array");
  }

  double gain = RecurseHelper(values.data(), static_cast<int>(values.size()),
                              cuts.data(), static_cast<int>(cuts.size()), which,
                              starts.data(), static_cast<int>(starts.size()),
                              resultCodes.data(), nPossibleRes);

  nb::list cutObj;
  for (auto cut : cuts) {
    cutObj.append(cut);
  }
  return nb::make_tuple(gain, cutObj);
}

static nb::list cQuantize_FindStartPoints(nb::object values, nb::object results,
                                          int nData) {
  nb::list startPts;
  if (nData < 2) {
    return startPts;
  }

  auto vals = sequenceToVector<double>(values);
  auto res = sequenceToVector<long int>(results);
  if (static_cast<size_t>(nData) > vals.size() ||
      static_cast<size_t>(nData) > res.size()) {
    throw nb::value_error("nData exceeds the length of values or results");
  }

  bool firstBlock = true;
  long lastBlockAct = -2, blockAct = res[0];
  int lastDiv = -1;
  constexpr double tol = 1e-8;

  int i = 1;
  while (i < nData) {
    while (i < nData && vals[i] - vals[i - 1] <= tol) {
      if (res[i] != blockAct) {
        blockAct = -1;
      }
      ++i;
    }
    if (firstBlock) {
      firstBlock = false;
      lastBlockAct = blockAct;
      lastDiv = i;
    } else {
      if (blockAct == -1 || lastBlockAct == -1 || blockAct != lastBlockAct) {
        startPts.append(lastDiv);
        lastDiv = i;
        lastBlockAct = blockAct;
      } else {
        lastDiv = i;
      }
    }
    if (i < nData) {
      blockAct = res[i];
    }
    ++i;
  }

  if (blockAct != lastBlockAct) {
    startPts.append(lastDiv);
  }
  return startPts;
}

NB_MODULE(cQuantize, m) {
  m.def("_RecurseOnBounds", cQuantize_RecurseOnBounds, "vals"_a, "pyCuts"_a,
        "which"_a, "pyStarts"_a, "results"_a, "nPossibleRes"_a,
        R"DOC(Find the best information gain for a set of quantization bounds.

Arguments:
  vals: sorted sequence of variable values
  pyCuts: indices into pyStarts
  which: index of the bound currently being adjusted
  pyStarts: potential starting points for quantization bounds
  results: sequence of integer result codes
  nPossibleRes: number of possible result codes

Returns:
  A tuple containing the best information gain and the corresponding cuts.
)DOC");
  m.def("_FindStartPoints", cQuantize_FindStartPoints, "values"_a, "results"_a,
        "nData"_a,
        R"DOC(Find possible quantization boundaries in sorted data.

Arguments:
  values: sorted sequence of variable values
  results: sequence of integer result codes
  nData: number of data points to examine

Returns:
  A list of indices that can be used as quantization boundaries.
)DOC");
}
