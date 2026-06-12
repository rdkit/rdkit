%include "std_vector.i"

%{
#include <vector>
#include <GraphMol/MolStandardize/Tautomer.h>
%}

%shared_ptr(RDKit::ROMol)
%template(Sizet_Vect) std::vector<size_t>;
%ignore RDKit::MolStandardize::TautomerScoringFunctions::makeOptimizedScorer;
%shared_ptr(RDKit::MolStandardize::TautomerEnumerator)
%newobject RDKit::MolStandardize::tautomerEnumeratorFromParams;
%newobject RDKit::MolStandardize::getV1TautomerEnumerator;
%newobject RDKit::MolStandardize::TautomerEnumerator::pickCanonical;
%newobject RDKit::MolStandardize::TautomerEnumerator::canonicalize;

%include <GraphMol/MolStandardize/Tautomer.h>

