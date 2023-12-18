
%include <std_unique_ptr.i>
%unique_ptr(ExplicitBitVect)

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

%extend RDKit::GeneralizedSubstruct::ExtendedQueryMol  {
  std::vector< std::vector<std::pair<int, int> > > getSubstructMatches(RDKit::ROMol &target,RDKit::SubstructMatchParameters ps = RDKit::SubstructMatchParameters()){
    std::vector<RDKit::MatchVectType> mvs = SubstructMatch(target, *($self),ps);
    return mvs;
  };
}
