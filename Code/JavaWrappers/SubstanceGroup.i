
%{
#include <GraphMol/SubstanceGroup.h>

RDKit::SubstanceGroup *getSubstanceGroupWithIdx(RDKit::ROMol &mol, unsigned int idx) {
  auto &groups = RDKit::getSubstanceGroups(mol);
  return &(groups[idx]);
}

unsigned int getSubstanceGroupCount(RDKit::ROMol &mol) {
    return RDKit::getSubstanceGroups(mol).size();
}
%}

// Base class RDProps is wrapped with shared_ptr, so SubstanceGroup must be too.
%shared_ptr(RDKit::SubstanceGroup)
%ignore getSubstanceGroups;

RDKit::SubstanceGroup *getSubstanceGroupWithIdx(RDKit::ROMol &mol, unsigned int idx);
unsigned int getSubstanceGroupCount(RDKit::ROMol &mol);

%include <GraphMol/SubstanceGroup.h>

%template(getStringProp) RDKit::SubstanceGroup::getProp<std::string>;
%template(getUIntProp) RDKit::SubstanceGroup::getProp<unsigned int>;
%template(getStringVectProp) RDKit::SubstanceGroup::getProp<RDKit::STR_VECT>;

