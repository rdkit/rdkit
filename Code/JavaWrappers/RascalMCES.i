
%include "std_vector.i"

%{
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
%}

#if SWIGCSHARP
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(Uint_Vect_Vect, std::vector<std::vector<unsigned int> >);
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(Uint_Vect, std::vector<unsigned int>);
#else
%template(Uint_Vect_Vect) std::vector<std::vector<unsigned int> >;
%template(Uint_Vect) std::vector<unsigned int>;
#endif
%template(Int_Vect_Vect) std::vector<std::vector<int> >;

%include <GraphMol/RascalMCES/RascalClusterOptions.h>
%include <GraphMol/RascalMCES/RascalOptions.h>
%include <GraphMol/RascalMCES/RascalResult.h>
%include <GraphMol/RascalMCES/RascalMCES.h>

