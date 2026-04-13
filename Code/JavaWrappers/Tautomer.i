%include "std_vector.i"

%{
#include <vector>
#include <GraphMol/MolStandardize/Tautomer.h>
%}

%shared_ptr(RDKit::ROMol)
%template(Sizet_Vect) std::vector<size_t>;
%ignore RDKit::MolStandardize::TautomerScoringFunctions::makeOptimizedScorer;

%include <GraphMol/MolStandardize/Tautomer.h>

