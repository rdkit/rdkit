//%import "ROMol.i"
%include "std_vector.i"


%{
#include <GraphMol/TautomerQuery/TautomerQuery.h>
%}
%shared_ptr(RDKit::TautomerQuery)
%shared_ptr(RDKit::ROMol)

%include <GraphMol/TautomerQuery/TautomerQuery.h>

