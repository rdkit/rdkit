//%import "ROMol.i"
%include "std_vector.i"


%{
#include <vector>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
%}
%shared_ptr(RDKit::ROMol)
%template(Sizet_Vect) std::vector<size_t>;

%include <GraphMol/TautomerQuery/TautomerQuery.h>

