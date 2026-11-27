#include <cstring>
#include <RDGeneral/Invariant.h>
#include <ML/InfoTheory/InfoGainFuncs.h>

/***********************************************

   constructs a variable table for the data passed in
   The table for a given variable records the number of times each possible
 value
    of that variable appears for each possible result of the function.

  **Arguments**

   - vals: pointer to double, contains the values of the variable,
     should be sorted

   - nVals: int, the length of _vals_

   - cuts: pointer to int, the indices of the quantization bounds

   - nCuts: int, the length of _cuts_

   - starts: pointer to int, the potential starting points for quantization
 bounds

   - nStarts: int, the length of _starts_

   - results: poitner to int, the result codes

   - nPossibleRes: int, the number of possible result codes


  **Returns**

    _varTable_ (a pointer to int), which is also modified in place.

  **Notes:**

    - _varTable_ is modified in place

    - the _results_ array is assumed to be _nVals_ long

 ***********************************************/
inline long int *GenVarTable(double *vals, int nVals, long int *cuts, int nCuts,
                             long int *starts, long int *results,
                             int nPossibleRes, long int *varTable) {
  RDUNUSED_PARAM(vals);
  int nBins = nCuts + 1;
  int idx, i, iTab;

  memset(varTable, 0, nBins * nPossibleRes * sizeof(long int));
  idx = 0;
  for (i = 0; i < nCuts; i++) {
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

/***********************************************

 This actually does the recursion required by *cQuantize_RecurseOnBounds()*,
  we do things this way to avoid having to convert things back and forth
  from Python objects

  **Arguments**

   - vals: pointer to double, contains the values of the variable,
     should be sorted

   - nVals: int, the length of _vals_

   - cuts: pointer to int, the indices of the quantization bounds

   - nCuts: int, the length of _cuts_

   - which: int, the quant bound being modified here

   - starts: pointer to int, the potential starting points for quantization
 bounds

   - nStarts: int, the length of _starts_

   - results: poitner to int, the result codes

   - nPossibleRes: int, the number of possible result codes


  **Returns**

    a double, the expected information gain for the best bounds found
      (which are found in _cuts_ )

  **Notes:**

    - _cuts_ is modified in place

    - the _results_ array is assumed to be _nVals_ long

 ***********************************************/
inline double RecurseHelper(double *vals, int nVals, long int *cuts, int nCuts,
                            int which, long int *starts, int nStarts,
                            long int *results, int nPossibleRes) {
  PRECONDITION(vals, "bad vals pointer");

  double maxGain = -1e6, gainHere;
  long int *bestCuts, *tCuts;
  long int *varTable = nullptr;
  int highestCutHere = nStarts - nCuts + which;
  int i, nBounds = nCuts;

  varTable = (long int *)calloc((nCuts + 1) * nPossibleRes, sizeof(long int));
  bestCuts = (long int *)calloc(nCuts, sizeof(long int));
  tCuts = (long int *)calloc(nCuts, sizeof(long int));
  CHECK_INVARIANT(varTable, "failed to allocate memory");
  CHECK_INVARIANT(bestCuts, "failed to allocate memory");
  CHECK_INVARIANT(tCuts, "failed to allocate memory");
  GenVarTable(vals, nVals, cuts, nCuts, starts, results, nPossibleRes,
              varTable);
  while (cuts[which] <= highestCutHere) {
    gainHere = RDInfoTheory::InfoEntropyGain(varTable, nCuts + 1, nPossibleRes);
    if (gainHere > maxGain) {
      maxGain = gainHere;
      memcpy(bestCuts, cuts, nCuts * sizeof(long int));
    }

    // recurse on the next vars if needed
    if (which < nBounds - 1) {
      memcpy(tCuts, cuts, nCuts * sizeof(long int));
      gainHere = RecurseHelper(vals, nVals, tCuts, nCuts, which + 1, starts,
                               nStarts, results, nPossibleRes);
      if (gainHere > maxGain) {
        maxGain = gainHere;
        memcpy(bestCuts, tCuts, nCuts * sizeof(long int));
      }
    }

    // update this cut
    int oldCut = cuts[which];
    cuts[which] += 1;
    int top, bot;
    bot = starts[oldCut];
    if (oldCut + 1 < nStarts) {
      top = starts[oldCut + 1];
    } else {
      top = starts[nStarts - 1];
    }
    for (i = bot; i < top; i++) {
      int v = results[i];
      varTable[which * nPossibleRes + v] += 1;
      varTable[(which + 1) * nPossibleRes + v] -= 1;
    }
    for (i = which + 1; i < nBounds; i++) {
      if (cuts[i] == cuts[i - 1]) {
        cuts[i] += 1;
      }
    }
  }
  memcpy(cuts, bestCuts, nCuts * sizeof(long int));
  free(tCuts);
  free(bestCuts);
  free(varTable);
  return maxGain;
}
