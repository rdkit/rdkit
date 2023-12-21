

#if SWIG_VERSION >= 0x040101                 
%include <std_unique_ptr.i>
%unique_ptr(ExplicitBitVect)
#endif

%{
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
%}
// %include "std_unique_ptr.i"
// %unique_ptr(ExtendedQueryMol)

#if SWIG_VERSION < 0x040101                 
%ignore patternFingerprintTargetMol(const ROMol &mol, unsigned int fpSize = 2048U);
%ignore patternFingerprintQuery(unsigned int fpSize = 2048U) const;
#endif

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
