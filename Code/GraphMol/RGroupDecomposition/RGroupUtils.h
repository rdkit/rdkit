#ifndef RGROUP_UTILS
#define RGROUP_UTILS

#include <GraphMol/RDKitBase.h>
#include <map>
namespace RDKit {
extern const std::string RLABEL;
extern const std::string SIDECHAIN_RLABELS;
extern const std::string done;

const unsigned int EMPTY_CORE_LABEL = -100000;

// Various places where rgroups can be labeled
//  the order of precedence
enum class Labelling {
  RGROUP_LABELS,
  ISOTOPE_LABELS,
  ATOMMAP_LABELS,
  INDEX_LABELS,
  DUMMY_LABELS,
  INTERNAL_LABELS
};

std::string labellingToString(Labelling type);
 
std::map<int, Atom *> getRlabels(const RWMol &mol);
void clearInputLabels(Atom *atom);
bool setLabel(Atom *atom, int label, std::set<int> &labels, int &maxLabel,
              bool relabel, Labelling type);
bool hasDummy(const RWMol &core);
}

#endif
