
%{
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
%}
// %include "std_unique_ptr.i"
// %unique_ptr(ExtendedQueryMol)

%ignore ExtendedQueryMol(std::unique_ptr<RWMol> mol);
%ignore ExtendedQueryMol(std::unique_ptr<MolBundle> mol);
%ignore ExtendedQueryMol(std::unique_ptr<TautomerQuery> mol);
%ignore ExtendedQueryMol(
      std::unique_ptr<std::vector<std::unique_ptr<TautomerQuery>>> tqs);
%ignore xqmol;

%include "GraphMol/GeneralizedSubstruct/XQMol.h";

